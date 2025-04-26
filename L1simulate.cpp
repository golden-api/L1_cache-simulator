#include <bits/stdc++.h>
using namespace std;

// Enums for MESI states, operation types, and bus operations
enum MESIState { INVALID,
                 SHARED,
                 EXCLUSIVE,
                 MODIFIED };
enum OpType { READ,
              WRITE };
enum BusOp { BUS_RD,
             BUS_RDX };

// Structure for a cache line
struct CacheLine {
  uint32_t tag;
  MESIState state;
  bool dirty;
  int lru_order;
  CacheLine() : tag(0), state(INVALID), dirty(false), lru_order(0) {}
};

// Structure for a bus transaction
struct BusTransaction {
  int core_id;
  BusOp op;
  uint32_t address;
  uint32_t block_addr;
  bool needs_writeback;  // Added to track if writeback is needed
  BusTransaction() : core_id(-1), op(BUS_RD), address(0), block_addr(0), needs_writeback(false) {}
  BusTransaction(int id, BusOp o, uint32_t addr, uint32_t baddr, bool wb = false)
      : core_id(id), op(o), address(addr), block_addr(baddr), needs_writeback(wb) {}
};

// Cache class to manage L1 cache
class Cache {
 private:
  int num_sets;
  int associativity;
  int block_size;
  int set_index_bits;
  int block_offset_bits;
  vector<vector<CacheLine>> cache_array;
  int miss_count;
  int access_count;
  int eviction_count;
  int writeback_count;

  void update_lru(int set_index, int way) {
    for (int i = 0; i < associativity; i++) {
      if (cache_array[set_index][i].state != INVALID) {
        cache_array[set_index][i].lru_order++;
      }
    }
    cache_array[set_index][way].lru_order = 0;
  }

  int find_lru_victim(int set_index) {
    int max_lru = -1;
    int victim = 0;
    for (int i = 0; i < associativity; i++) {
      if (cache_array[set_index][i].state == INVALID) return i;
      if (cache_array[set_index][i].lru_order > max_lru) {
        max_lru = cache_array[set_index][i].lru_order;
        victim = i;
      }
    }
    return victim;
  }

 public:
  Cache(int s, int E, int b) : set_index_bits(s), associativity(E), block_offset_bits(b) {
    num_sets = 1 << s;
    block_size = 1 << b;
    cache_array.resize(num_sets, vector<CacheLine>(associativity));
    miss_count = access_count = eviction_count = writeback_count = 0;
  }

  bool access(uint32_t address, OpType op, bool& needs_bus, BusOp& bus_op, uint32_t& block_addr, bool& needs_writeback) {
    access_count++;
    uint32_t set_index = (address >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = address >> (block_offset_bits + set_index_bits);
    block_addr = address & ~((1 << block_offset_bits) - 1);

    for (int i = 0; i < associativity; i++) {
      CacheLine& line = cache_array[set_index][i];
      if (line.state != INVALID && line.tag == tag) {
        update_lru(set_index, i);
        if (op == READ) {
          needs_bus = false;
          return true;
        } else {
          if (line.state == MODIFIED || line.state == EXCLUSIVE) {
            line.dirty = true;
            line.state = MODIFIED;
            needs_bus = false;
            return true;
          } else if (line.state == SHARED) {
            needs_bus = true;
            bus_op = BUS_RDX;
            return true;
          }
        }
      }
    }

    miss_count++;
    needs_bus = true;
    bus_op = (op == READ) ? BUS_RD : BUS_RDX;

    int victim = find_lru_victim(set_index);
    CacheLine& line = cache_array[set_index][victim];
    if (line.state != INVALID) {
      eviction_count++;
      if (line.state == MODIFIED) {
        writeback_count++;
        needs_writeback = true;  // Set writeback flag for the transaction
      }
    }

    line.tag = tag;
    line.state = (op == READ) ? EXCLUSIVE : MODIFIED;
    line.dirty = (op == WRITE);
    update_lru(set_index, victim);
    return false;
  }

  bool snoop(const BusTransaction& trans) {  // Returns true if invalidation occurred
    uint32_t set_index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = trans.block_addr >> (block_offset_bits + set_index_bits);
    for (auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) {
        if (trans.op == BUS_RD) {
          if (line.state == EXCLUSIVE || line.state == MODIFIED) {
            line.state = SHARED;
          }
        } else if (trans.op == BUS_RDX) {
          line.state = INVALID;
          return true;  // Invalidation occurred
        }
        break;
      }
    }
    return false;
  }

  void update_after_bus(const BusTransaction& trans) {
    uint32_t set_index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = trans.block_addr >> (block_offset_bits + set_index_bits);
    bool found = false;
    for (int i = 0; i < associativity; i++) {
      CacheLine& line = cache_array[set_index][i];
      if (line.state != INVALID && line.tag == tag) {
        found = true;
        if (trans.op == BUS_RDX) {
          line.state = MODIFIED;
          line.dirty = true;
        }
        update_lru(set_index, i);
        break;
      }
    }
    if (!found && (trans.op == BUS_RD || trans.op == BUS_RDX)) {
      int way = -1;
      for (int i = 0; i < associativity; i++) {
        if (cache_array[set_index][i].state == INVALID) {
          way = i;
          break;
        }
      }
      if (way == -1) {
        way = find_lru_victim(set_index);
        CacheLine& victim = cache_array[set_index][way];
        if (victim.state != INVALID) {
          eviction_count++;
          if (victim.state == MODIFIED) {
            writeback_count++;
          }
        }
      }
      CacheLine& line = cache_array[set_index][way];
      line.tag = tag;
      line.state = (trans.op == BUS_RD) ? SHARED : MODIFIED;
      line.dirty = (trans.op == BUS_RDX);
      update_lru(set_index, way);
    }
  }

  bool has_block(uint32_t block_addr) const {
    uint32_t set_index = (block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = block_addr >> (block_offset_bits + set_index_bits);
    for (const auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) {
        return true;
      }
    }
    return false;
  }

  int get_miss_count() const { return miss_count; }
  int get_access_count() const { return access_count; }
  int get_eviction_count() const { return eviction_count; }
  int get_writeback_count() const { return writeback_count; }
  int get_block_size() const { return block_size; }
};

// Bus class to manage bus transactions
class Bus {
 private:
  BusTransaction current_transaction;
  int remaining_cycles;
  queue<BusTransaction> pending;
  int invalidations;
  uint64_t data_traffic;

 public:
  Bus() : remaining_cycles(0), invalidations(0), data_traffic(0), current_transaction() {}

  void add_transaction(const BusTransaction& trans) {
    pending.push(trans);
  }

  void start_next_transaction(const vector<Cache*>& caches) {
    if (!pending.empty()) {
      current_transaction = pending.front();
      pending.pop();
      uint32_t block_addr = current_transaction.block_addr;
      bool found = false;
      for (int i = 0; i < caches.size(); i++) {
        if (i == current_transaction.core_id) continue;
        if (caches[i]->has_block(block_addr)) {
          found = true;
          break;
        }
      }
      int block_size = caches[0]->get_block_size();
      int N = block_size / 4;
      if (found) {
        remaining_cycles = 2 * N;
      } else {
        remaining_cycles = 100;
      }
      if (current_transaction.needs_writeback) {
        remaining_cycles += 100;  // Add writeback time
      }
      data_traffic += block_size;
    } else {
      remaining_cycles = 0;
    }
  }

  bool process_cycle() {
    if (remaining_cycles > 0) {
      remaining_cycles--;
      if (remaining_cycles == 0) {
        return true;
      }
    }
    return false;
  }

  bool is_busy() const { return remaining_cycles > 0; }

  const BusTransaction& get_current_transaction() const { return current_transaction; }

  int get_invalidations() const { return invalidations; }

  uint64_t get_data_traffic() const { return data_traffic; }

  void increment_invalidations() { invalidations++; }

  bool is_pending_empty() const { return pending.empty(); }

  bool has_pending() const { return !pending.empty(); }
};

// Core class to manage a processor core
class Core {
 private:
  int id;
  ifstream trace_file;
  int instruction_count;
  uint64_t total_cycles;
  uint64_t idle_cycles;
  Cache* cache;
  bool is_stalled;

 public:
  Core(int id, const string& trace_filename, int s, int E, int b)
      : id(id), instruction_count(0), total_cycles(0), idle_cycles(0), is_stalled(false) {
    cache = new Cache(s, E, b);
    trace_file.open(trace_filename);
    if (!trace_file.is_open()) {
      cerr << "Error opening trace file: " << trace_filename << endl;
      exit(1);
    }
  }

  ~Core() {
    delete cache;
    if (trace_file.is_open()) trace_file.close();
  }

  bool has_next_instruction() {
    return trace_file.peek() != EOF;
  }

  void process_instruction(Bus* bus) {
    string line;
    getline(trace_file, line);
    istringstream iss(line);
    char op_char;
    uint32_t addr;
    iss >> op_char >> hex >> addr;
    OpType op_type = (op_char == 'R') ? READ : WRITE;
    bool needs_bus;
    BusOp bus_op;
    uint32_t block_addr;
    bool needs_writeback = false;
    bool hit = cache->access(addr, op_type, needs_bus, bus_op, block_addr, needs_writeback);
    if (!hit || needs_bus) {
      BusTransaction trans(id, bus_op, addr, block_addr, needs_writeback);
      bus->add_transaction(trans);
      is_stalled = true;
    }
    instruction_count++;
  }

  void set_stalled(bool stalled) { is_stalled = stalled; }

  void increment_idle_cycles() { idle_cycles++; }

  void set_total_cycles(uint64_t cycles) { total_cycles = cycles; }

  int get_id() const { return id; }
  int get_instruction_count() const { return instruction_count; }
  uint64_t get_total_cycles() const { return total_cycles; }
  uint64_t get_idle_cycles() const { return idle_cycles; }
  double get_miss_rate() const {
    return cache->get_access_count() > 0 ? static_cast<double>(cache->get_miss_count()) / cache->get_access_count() : 0.0;
  }
  int get_eviction_count() const { return cache->get_eviction_count(); }
  int get_writeback_count() const { return cache->get_writeback_count(); }
  Cache* get_cache() { return cache; }
  bool get_is_stalled() const { return is_stalled; }
};

// Simulator class to orchestrate the simulation
class Simulator {
 private:
  vector<Core*> cores;
  Bus* bus;
  string app_name;
  int set_index_bits;
  int assoc;
  int block_offset_bits;  // Renamed from block_bits for clarity
  string output_file;
  uint64_t global_cycle;

  vector<Cache*> get_caches() {
    vector<Cache*> caches;
    for (auto& core : cores) caches.push_back(core->get_cache());
    return caches;
  }

 public:
  Simulator(const string& app, int s, int E, int b, const string& out)
      : app_name(app), set_index_bits(s), assoc(E), block_offset_bits(b), output_file(out), global_cycle(0) {
    bus = new Bus();
    for (int i = 0; i < 4; i++) {
      string trace = app_name + "_proc" + to_string(i) + ".trace";
      cores.push_back(new Core(i, trace, s, E, b));
    }
  }

  ~Simulator() {
    for (Core* core : cores) delete core;
    delete bus;
  }

  void run() {
    while (true) {
      bool all_done = true;
      for (auto& core : cores) {
        if (core->has_next_instruction() || core->get_is_stalled()) {
          all_done = false;
          break;
        }
      }
      if (all_done && !bus->is_busy() && bus->is_pending_empty()) break;

      // Process bus transaction for one cycle
      if (bus->process_cycle()) {
        int core_id = bus->get_current_transaction().core_id;
        if (core_id >= 0 && core_id < 4) {
          cores[core_id]->set_stalled(false);
          cores[core_id]->get_cache()->update_after_bus(bus->get_current_transaction());
        }
        for (int i = 0; i < 4; i++) {
          if (i != core_id) {
            Cache* cache = cores[i]->get_cache();
            if (cache->snoop(bus->get_current_transaction())) {
              bus->increment_invalidations();
            }
          }
        }
      }

      // Start next transaction if bus is idle
      if (!bus->is_busy() && bus->has_pending()) {
        bus->start_next_transaction(get_caches());
      }

      // Process each core
      for (auto& core : cores) {
        if (!core->get_is_stalled() && core->has_next_instruction()) {
          core->process_instruction(bus);
        }
        if (core->get_is_stalled()) {
          core->increment_idle_cycles();
        }
      }

      global_cycle++;
    }

    // Set total cycles for each core
    for (auto& core : cores) {
      core->set_total_cycles(global_cycle);
    }
  }

  void write_output() {
    ofstream out(output_file);
    for (int i = 0; i < 4; i++) {
      out << "Core " << i << ":\n";
      out << "Instructions: " << cores[i]->get_instruction_count() << "\n";
      out << "Total Cycles: " << cores[i]->get_total_cycles() << "\n";
      out << "Idle Cycles: " << cores[i]->get_idle_cycles() << "\n";
      out << "Miss Rate: " << fixed << setprecision(4) << cores[i]->get_miss_rate() << "\n";
      out << "Evictions: " << cores[i]->get_eviction_count() << "\n";
      out << "Writebacks: " << cores[i]->get_writeback_count() << "\n";
    }
    out << "System:\n";
    out << "Invalidations: " << bus->get_invalidations() << "\n";
    out << "Bus Traffic: " << bus->get_data_traffic() << " bytes\n";
    out.close();
  }
};

// Command-line argument parser
struct Args {
  string app_name;
  int set_index_bits;
  int assoc;
  int block_offset_bits;  // Renamed for clarity
  string output_file;
  bool show_help;
};

void print_help() {
  cout << "./L1simulate\n"
       << "-t <tracefile>: name of parallel application (e.g., app1)\n"
       << "-s <s>: number of set index bits\n"
       << "-E <E>: associativity\n"
       << "-b <b>: number of block offset bits\n"
       << "-o <outfilename>: output file\n"
       << "-h: prints this help\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args = {"", 0, 0, 0, "", false};
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if (arg == "-h") {
      args.show_help = true;
      return args;
    }
    if (i + 1 >= argc) {
      cerr << "Missing value for option " << arg << endl;
      args.show_help = true;
      return args;
    }
    string value = argv[i + 1];
    if (arg == "-t")
      args.app_name = value;
    else if (arg == "-s")
      args.set_index_bits = stoi(value);
    else if (arg == "-E")
      args.assoc = stoi(value);
    else if (arg == "-b")
      args.block_offset_bits = stoi(value);  // Updated
    else if (arg == "-o")
      args.output_file = value;
    else {
      cerr << "Unknown option: " << arg << endl;
      args.show_help = true;
    }
    ++i;
  }
  return args;
}

int main(int argc, char* argv[]) {
  Args args = parse_args(argc, argv);
  if (args.show_help || args.app_name.empty() || args.set_index_bits <= 0 ||
      args.assoc <= 0 || args.block_offset_bits <= 0 || args.output_file.empty()) {
    print_help();
    return 1;
  }
  Simulator sim(args.app_name, args.set_index_bits, args.assoc, args.block_offset_bits, args.output_file);
  sim.run();
  sim.write_output();
  return 0;
}