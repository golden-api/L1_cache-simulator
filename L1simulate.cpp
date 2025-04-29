#include <bits/stdc++.h>
using namespace std;
#define int uint32_t

// Enums for clarity in MESI states, operation types, and bus operations
enum MESIState { INVALID,
                 SHARED,
                 EXCLUSIVE,
                 MODIFIED };
enum OpType { READ,
              WRITE };
enum BusOp { BUS_RD,
             BUS_RDX };

// CacheLine Structure: Represents a single line in the cache
struct CacheLine {
  int tag;
  MESIState state;
  bool dirty;
  int lru_order;
  CacheLine() : tag(0), state(INVALID), dirty(false), lru_order(0) {}
};

// BusTransaction Structure: Represents a transaction on the bus
struct BusTransaction {
  int core_id;      // Requestor core ID
  int receptor_id;  // Receptor core ID (core with Modified line), -1 if none
  BusOp op;
  int address;
  int block_addr;
  bool needs_writeback;
  bool found_in_other;
  BusTransaction() : core_id(-1), receptor_id(-1), op(BUS_RD), address(0), block_addr(0), needs_writeback(false), found_in_other(false) {}
  BusTransaction(int id, BusOp o, int addr, int baddr, bool wb = false)
      : core_id(id), receptor_id(-1), op(o), address(addr), block_addr(baddr), needs_writeback(wb), found_in_other(false) {}
};

// Forward declaration of Core class
class Core;

// Cache Class: Manages an L1 cache for one core
class Cache {
 private:
  int num_sets;
  int associativity;
  int block_size;
  int set_index_bits;
  int block_offset_bits;
  vector<vector<CacheLine>> cache_array;
  int access_count;
  int miss_count;
  int eviction_count;
  int writeback_count;
  int invalidation_count;  // Tracks coherence invalidations

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
    access_count = miss_count = eviction_count = writeback_count = invalidation_count = 0;
  }

  bool access(int address, OpType op, bool& needs_bus, BusOp& bus_op, int& block_addr, bool& needs_writeback) {
    access_count++;

    int set_index = (address >> block_offset_bits) & ((1u << set_index_bits) - 1);
    int tag = address >> (block_offset_bits + set_index_bits);
    block_addr = address & ~((1u << block_offset_bits) - 1);

    for (int way = 0; way < associativity; way++) {
      CacheLine& line = cache_array[set_index][way];
      if (line.state != INVALID && line.tag == tag) {
        update_lru(set_index, way);
        if (op == READ) {
          needs_bus = false;
          return true;
        }
        switch (line.state) {
          case MODIFIED:
            line.dirty = true;
            needs_bus = false;
            return true;
          case EXCLUSIVE:
            line.state = MODIFIED;
            line.dirty = true;
            needs_bus = false;
            return true;
          case SHARED:
            needs_bus = true;
            bus_op = BUS_RDX;
            return true;
          default:
            break;
        }
      }
    }

    miss_count++;
    needs_bus = true;
    bus_op = (op == READ ? BUS_RD : BUS_RDX);
    needs_writeback = false;

    int victim = find_lru_victim(set_index);
    CacheLine& victim_line = cache_array[set_index][victim];
    if (victim_line.state != INVALID) {
      eviction_count++;
      if (victim_line.state == MODIFIED) {
        writeback_count++;
        needs_writeback = true;
      }
    }

    return false;
  }

  bool snoop(const BusTransaction& trans) {
    int set_index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    int tag = trans.block_addr >> (block_offset_bits + set_index_bits);

    for (auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) {
        if (trans.op == BUS_RD) {
          if (line.state == EXCLUSIVE || line.state == MODIFIED) {
            if (line.state == MODIFIED) {
              writeback_count++;
            }
            line.state = SHARED;
          }
        } else if (trans.op == BUS_RDX) {
          if (line.state != INVALID) {
            if (line.state == MODIFIED) {
              writeback_count++;
            }
            line.state = INVALID;
            invalidation_count++;
            return true;
          }
        }
        break;
      }
    }
    return false;
  }

  void update_after_bus(const BusTransaction& trans) {
    int set_index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    int tag = trans.block_addr >> (block_offset_bits + set_index_bits);
    bool found = false;
    for (int i = 0; i < associativity; i++) {
      CacheLine& line = cache_array[set_index][i];
      if (line.state != INVALID && line.tag == tag) {
        found = true;
        if (trans.op == BUS_RD) {
          line.state = trans.found_in_other ? SHARED : EXCLUSIVE;
        } else if (trans.op == BUS_RDX) {
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
      if (trans.op == BUS_RD) {
        line.state = trans.found_in_other ? SHARED : EXCLUSIVE;
      } else if (trans.op == BUS_RDX) {
        line.state = MODIFIED;
        line.dirty = true;
      }
      update_lru(set_index, way);
    }
  }

  bool has_block(int block_addr) const {
    int set_index = (block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    int tag = block_addr >> (block_offset_bits + set_index_bits);
    for (const auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) return true;
    }
    return false;
  }

  MESIState get_line_state(int block_addr) const {
    int set_index = (block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    int tag = block_addr >> (block_offset_bits + set_index_bits);
    for (const auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) {
        return line.state;
      }
    }
    return INVALID;
  }

  int get_miss_count() const { return miss_count; }
  int get_access_count() const { return access_count; }
  int get_eviction_count() const { return eviction_count; }
  int get_writeback_count() const { return writeback_count; }
  int get_invalidation_count() const { return invalidation_count; }
  int get_block_size() const { return block_size; }
};

// Bus Class: Manages bus transactions with separate cycle attribution
class Bus {
 private:
  BusTransaction current_transaction;
  int remaining_cycles;
  deque<BusTransaction> pending;
  vector<BusTransaction> cycle_transactions;
  int invalidations;
  uint64_t data_traffic;
  mt19937 rng;

 public:
  int total_transactions = 0;
  unordered_set<int> pending_read_blocks;
  unordered_map<int, vector<int>> rd_waiters;

  Bus() : remaining_cycles(0), invalidations(0), data_traffic(0), current_transaction() {
    random_device rd;
    rng = mt19937(rd());
  }

  void add_transaction(const BusTransaction& t) {
    ++total_transactions;
    if (t.op == BUS_RD) {
      auto block = t.block_addr;
      if (pending_read_blocks.count(block)) {
        rd_waiters[block].push_back(t.core_id);
        return;
      }
      pending_read_blocks.insert(block);
    }
    cycle_transactions.push_back(t);
  }

  void commit_transactions() {
    for (auto& t : cycle_transactions) {
      pending.push_back(t);
    }
    cycle_transactions.clear();
  }

  void start_next_transaction(const vector<Cache*>& caches) {
    if (pending.empty()) return;
    current_transaction = pending.front();
    pending.pop_front();

    int block = current_transaction.block_addr;
    bool found = false;
    bool foundM = false;
    int receptor_id = -1;
    for (int i = 0; i < (int)caches.size(); i++) {
      if (i == current_transaction.core_id) continue;
      if (caches[i]->has_block(block)) {
        found = true;
        MESIState st = caches[i]->get_line_state(block);
        if (st == MODIFIED) {
          foundM = true;
          receptor_id = i;
          break;
        }
      }
    }
    current_transaction.found_in_other = found;
    current_transaction.receptor_id = receptor_id;

    int block_size = caches[0]->get_block_size();
    int N = block_size / 4;

    if (current_transaction.op == BUS_RD) {
      if (!found) {
        remaining_cycles = 100;
      } else {
        remaining_cycles = (foundM ? 100 : 0) + 2 * N;
      }
    } else if (current_transaction.op == BUS_RDX) {
      if (found) {
        if (foundM) {
          remaining_cycles = 201;  // 100 for writeback (receptor) + 101 for memory read (requestor)
        } else {
          remaining_cycles = 101;
        }
      } else {
        remaining_cycles = 101;
      }
      if (current_transaction.needs_writeback)
        remaining_cycles += 100;
    }

    data_traffic += block_size;
  }

  // Declare process_cycle; implementation moved after Core
  bool process_cycle(vector<Core*>& cores);

  bool is_busy() const { return remaining_cycles > 0; }
  BusTransaction& get_current_transaction() { return current_transaction; }
  int get_invalidations() const { return invalidations; }
  uint64_t get_data_traffic() const { return data_traffic; }
  void increment_invalidations() { invalidations++; }
  bool is_pending_empty() const { return pending.empty(); }
  bool has_pending() const { return !pending.empty(); }
};

// Core Class: Represents a processor core
class Core {
 private:
  int id;
  ifstream trace_file;
  int instruction_count;
  uint64_t total_cycles;
  uint64_t idle_cycles;
  Cache* cache;
  bool is_stalled;
  int read_count = 0, write_count = 0;

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

  bool has_next_instruction() { return trace_file.peek() != EOF; }

  void process_instruction(Bus* bus) {
    string line;
    getline(trace_file, line);
    istringstream iss(line);
    char op_char;
    int addr;
    iss >> op_char >> hex >> addr;
    OpType op_type = (op_char == 'R') ? READ : WRITE;
    bool needs_bus;
    BusOp bus_op;
    int block_addr;
    bool needs_writeback = false;
    if (op_type == READ) {
      ++read_count;
    } else {
      ++write_count;
    }
    bool hit = cache->access(addr, op_type, needs_bus, bus_op, block_addr, needs_writeback);
    total_cycles++;
    if (hit) {
      cout << "Core " << id << ": Cache hit for address 0x" << hex << addr << dec << endl;
    }
    if (!hit || needs_bus) {
      BusTransaction trans(id, bus_op, addr, block_addr, needs_writeback);
      bus->add_transaction(trans);
      is_stalled = true;
    }
    instruction_count++;
  }

  void set_stalled(bool stalled) { is_stalled = stalled; }
  void increment_idle_cycles() { idle_cycles++, total_cycles++; }
  void increment_total_cycles() { total_cycles++; }
  int get_id() const { return id; }
  int get_instruction_count() const { return instruction_count; }
  uint64_t get_total_cycles() const { return total_cycles; }
  uint64_t get_idle_cycles() const { return idle_cycles; }
  double get_miss_rate() const {
    return cache->get_access_count() > 0 ? static_cast<double>(cache->get_miss_count()) / cache->get_access_count() : 0.0;
  }
  int get_read_count() const { return read_count; }
  int get_write_count() const { return write_count; }
  int get_eviction_count() const { return cache->get_eviction_count(); }
  int get_writeback_count() const { return cache->get_writeback_count(); }
  int get_invalidation_count() const { return cache->get_invalidation_count(); }
  Cache* get_cache() { return cache; }
  bool get_is_stalled() const { return is_stalled; }
};

// Implementation of Bus::process_cycle (moved after Core)
bool Bus::process_cycle(vector<Core*>& cores) {
  if (remaining_cycles > 0) {
    remaining_cycles--;
    if (remaining_cycles == 0) {
      cout << "Bus transaction completed: core " << current_transaction.core_id
           << ", op " << (current_transaction.op == BUS_RD ? "BUS_RD" : "BUS_RDX")
           << ", address 0x" << hex << current_transaction.address << dec << endl;
    }
    return remaining_cycles == 0;
  }
  return false;
}

// Simulator Class: Orchestrates the simulation
class Simulator {
 private:
  vector<Core*> cores;
  Bus* bus;
  string app_name;
  int set_index_bits;
  int assoc;
  int block_offset_bits;
  string output_file;

  vector<Cache*> get_caches() {
    vector<Cache*> caches;
    for (auto& core : cores) caches.push_back(core->get_cache());
    return caches;
  }

 public:
  Simulator(const string& app, int s, int E, int b, const string& out)
      : app_name(app), set_index_bits(s), assoc(E), block_offset_bits(b), output_file(out) {
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
    int start_core = 0;
    while (true) {
      if (bus->is_busy()) {
        if (bus->process_cycle(cores)) {
          auto& done = bus->get_current_transaction();
          auto block = done.block_addr;

          // Update caches and unstall requesting core
          cores[done.core_id]->set_stalled(false);
          cores[done.core_id]->get_cache()->update_after_bus(done);

          if (done.op == BUS_RD) {
            auto it = bus->rd_waiters.find(block);
            if (it != bus->rd_waiters.end()) {
              for (int cid : it->second) {
                cores[cid]->set_stalled(false);
                cores[cid]->get_cache()->update_after_bus(done);
              }
              bus->rd_waiters.erase(it);
            }
          }

          for (int i = 0; i < 4; i++) {
            if (i != done.core_id && cores[i]->get_cache()->snoop(done)) {
              bus->increment_invalidations();
            }
          }

          if (done.op == BUS_RD) {
            bus->pending_read_blocks.erase(block);
          }
        }
      }

      if (!bus->is_busy() && bus->has_pending()) {
        bus->start_next_transaction(get_caches());
      }

      for (int i = 0; i < 4; i++) {
        int cid = (start_core + i) % 4;
        Core* core = cores[cid];

        if (!core->get_is_stalled() && core->has_next_instruction()) {
          core->process_instruction(bus);
        } else if (core->get_is_stalled()) {
          if (bus->is_busy() && cid == bus->get_current_transaction().core_id) {
            core->increment_total_cycles();
          } else {
            core->increment_idle_cycles();
          }
        }
      }

      bus->commit_transactions();
      start_core = (start_core + 1) % 4;

      bool all_done = true;
      for (auto& core : cores) {
        if (core->has_next_instruction() || core->get_is_stalled()) {
          all_done = false;
          break;
        }
      }
      if (all_done && !bus->is_busy() && bus->is_pending_empty()) {
        break;
      }
    }
  }
  // in class Simulator
  void write_output() {
    ofstream out(output_file);

    // Simulation parameters
    out << "Simulation Parameters:\n";
    out << "Trace Prefix: " << app_name << "\n";
    out << "Set Index Bits: " << set_index_bits << "\n";
    out << "Associativity: " << assoc << "\n";
    out << "Block Bits: " << block_offset_bits << "\n";
    out << "Block Size (Bytes): " << (1u << block_offset_bits) << "\n";
    out << "Number of Sets: " << (1u << set_index_bits) << "\n";
    out << "Cache Size (KB per core): "
        << ((1u << block_offset_bits) * assoc * (1u << set_index_bits) / 1024) << "\n";
    out << "MESI Protocol: Enabled\n";
    out << "Write Policy: Write-back, Write-allocate\n";
    out << "Replacement Policy: LRU\n";
    out << "Bus: Central snooping bus\n\n";

    // Per-core statistics
    for (int i = 0; i < (int)cores.size(); i++) {
      Core* c = cores[i];
      out << "Core " << i << " Statistics:\n";
      out << "Total Instructions: " << c->get_instruction_count() << "\n";
      out << "Total Reads: " << c->get_read_count() << "\n";
      out << "Total Writes: " << c->get_write_count() << "\n";
      out << "Total Execution Cycles: " << c->get_total_cycles() << "\n";
      out << "Idle Cycles: " << c->get_idle_cycles() << "\n";
      out << "Cache Misses: " << c->get_cache()->get_miss_count() << "\n";
      out << "Cache Miss Rate: " << fixed << setprecision(2)
          << (100.0 * c->get_cache()->get_miss_count() / c->get_cache()->get_access_count())
          << "%\n";
      out << "Cache Evictions: " << c->get_cache()->get_eviction_count() << "\n";
      out << "Writebacks: " << c->get_cache()->get_writeback_count() << "\n";
      out << "Bus Invalidations: " << c->get_cache()->get_invalidation_count() << "\n";
      out << "Data Traffic (Bytes): " << c->get_cache()->get_block_size() * /*assume*/ 1 /*#blocks transferred*/
          << "\n\n";
    }

    // Overall bus summary
    out << "Overall Bus Summary:\n";
    out << "Total Bus Transactions: " << bus->total_transactions << "\n";
    out << "Total Bus Traffic (Bytes): " << bus->get_data_traffic() << "\n";

    out.close();
  }
};

// Command-line argument parsing
struct Args {
  string app_name;
  int set_index_bits;
  int assoc;
  int block_offset_bits;
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
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-h") {
      args.show_help = true;
      return args;
    }
    if (i + 1 >= argc) {
      cerr << "Missing value for " << arg << endl;
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
      args.block_offset_bits = stoi(value);
    else if (arg == "-o")
      args.output_file = value;
    else {
      cerr << "Unknown option: " << arg << endl;
      args.show_help = true;
    }
    i++;
  }
  return args;
}

signed main(signed argc, char* argv[]) {
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