#include <bits/stdc++.h>
using namespace std;

// **Enums for clarity in MESI states, operation types, and bus operations**
enum MESIState { INVALID,
                 SHARED,
                 EXCLUSIVE,
                 MODIFIED };  // MESI protocol states
enum OpType { READ,
              WRITE };  // Types of memory operations
enum BusOp { BUS_RD,
             BUS_RDX };  // Bus operations: read (BusRd), read-exclusive (BusRdX)

// **CacheLine Structure: Represents a single line in the cache**
struct CacheLine {
  uint32_t tag;     // Tag bits to identify the memory block
  MESIState state;  // Current MESI state (INVALID, SHARED, EXCLUSIVE, MODIFIED)
  bool dirty;       // Dirty bit for write-back policy
  int lru_order;    // LRU order for replacement (0 = most recently used)
  CacheLine() : tag(0), state(INVALID), dirty(false), lru_order(0) {}
};

// **BusTransaction Structure: Represents a transaction on the bus**
struct BusTransaction {
  int core_id;           // ID of the core issuing the transaction
  BusOp op;              // Type of bus operation (BUS_RD or BUS_RDX)
  uint32_t address;      // Original memory address
  uint32_t block_addr;   // Block-aligned address
  bool needs_writeback;  // True if evicting a dirty block
  bool found_in_other;   // True if block was found in another cache
  BusTransaction() : core_id(-1), op(BUS_RD), address(0), block_addr(0), needs_writeback(false), found_in_other(false) {}
  BusTransaction(int id, BusOp o, uint32_t addr, uint32_t baddr, bool wb = false)
      : core_id(id), op(o), address(addr), block_addr(baddr), needs_writeback(wb), found_in_other(false) {}
};

// **Cache Class: Manages an L1 cache for one core**
class Cache {
 private:
  int num_sets;                           // Number of sets (S = 2^s)
  int associativity;                      // Number of ways (E)
  int block_size;                         // Block size in bytes (B = 2^b)
  int set_index_bits;                     // Number of set index bits (s)
  int block_offset_bits;                  // Number of block offset bits (b)
  vector<vector<CacheLine>> cache_array;  // 2D array: num_sets x associativity
  int miss_count;                         // Number of cache misses
  int access_count;                       // Total number of cache accesses
  int eviction_count;                     // Number of evictions
  int writeback_count;                    // Number of writebacks to memory

  // **Update LRU order: Set accessed line to MRU (0)**
  void update_lru(int set_index, int way) {
    for (int i = 0; i < associativity; i++) {
      if (cache_array[set_index][i].state != INVALID) {
        cache_array[set_index][i].lru_order++;
      }
    }
    cache_array[set_index][way].lru_order = 0;
  }

  // **Find LRU victim for replacement**
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

  // **Access cache: Handle read/write operations with independent case handling**
  bool access(uint32_t address, OpType op, bool& needs_bus, BusOp& bus_op, uint32_t& block_addr, bool& needs_writeback) {
    access_count++;

    // Decode index & tag
    uint32_t set_index = (address >> block_offset_bits) & ((1u << set_index_bits) - 1);
    uint32_t tag = address >> (block_offset_bits + set_index_bits);
    block_addr = address & ~((1u << block_offset_bits) - 1);

    // 1) Search for a hit
    for (int way = 0; way < associativity; way++) {
      CacheLine& line = cache_array[set_index][way];
      if (line.state != INVALID && line.tag == tag) {
        // Hit — update LRU
        update_lru(set_index, way);

        if (op == READ) {
          // --- Read hit: nothing on bus, no state change
          needs_bus = false;
          return true;
        }

        // --- WRITE hit: dispatch based on current MESI state
        switch (line.state) {
          case MODIFIED:
            // (M) already exclusive & dirty
            // Just update the value locally
            line.dirty = true;
            needs_bus = false;
            return true;

          case EXCLUSIVE:
            // (E) clean exclusive
            // Upgrade to Modified, mark dirty
            line.state = MODIFIED;
            line.dirty = true;
            needs_bus = false;
            return true;

          case SHARED:
            // (S) must invalidate all other copies
            // Broadcast BusRdX (RWITM) → other caches will snoop and go I
            needs_bus = true;
            bus_op = BUS_RDX;
            // leave state change & dirty bit to update_after_bus()
            return true;

          default:
            // INVALID should never match here
            break;
        }
      }
    }

    // 2) MISS handling (read or write) stays the same as before...
    //    (install deferred until update_after_bus)
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

    // **Do not install a new line here** — update_after_bus() will allocate.
    return false;
  }

  // **Snoop: Handle bus transactions from other caches**
  bool snoop(const BusTransaction& trans) {
    uint32_t set_index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = trans.block_addr >> (block_offset_bits + set_index_bits);

    for (auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) {
        if (trans.op == BUS_RD) {
          // Any E or M downgrades to SHARED
          if (line.state == EXCLUSIVE || line.state == MODIFIED) {
            if (line.state == MODIFIED) {
              // Write back the dirty data before sharing
              writeback_count++;
            }
            line.state = SHARED;
          }
        } else if (trans.op == BUS_RDX) {
          // Invalidate on a read-with-intent-to-modify
          if (line.state == MODIFIED) {
            writeback_count++;
          }
          line.state = INVALID;
          return true;
        }
        break;
      }
    }
    return false;
  }

  // **Update cache state after bus transaction**
  void update_after_bus(const BusTransaction& trans) {
    uint32_t set_index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = trans.block_addr >> (block_offset_bits + set_index_bits);
    bool found = false;
    for (int i = 0; i < associativity; i++) {
      CacheLine& line = cache_array[set_index][i];
      if (line.state != INVALID && line.tag == tag) {
        found = true;
        if (trans.op == BUS_RD) {
          line.state = trans.found_in_other ? SHARED : EXCLUSIVE;
          // cout << "Update after BUS_RD at block 0x" << hex << trans.block_addr << dec
          //      << ": state set to " << (trans.found_in_other ? "SHARED" : "EXCLUSIVE") << endl;
        } else if (trans.op == BUS_RDX) {
          line.state = MODIFIED;
          line.dirty = true;
          // cout << "Update after BUS_RDX at block 0x" << hex << trans.block_addr << dec
          //      << ": state set to MODIFIED" << endl;
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
            // cout << "Evicting MODIFIED line during update_after_bus, writeback required" << endl;
          }
        }
      }
      CacheLine& line = cache_array[set_index][way];
      line.tag = tag;
      if (trans.op == BUS_RD) {
        line.state = trans.found_in_other ? SHARED : EXCLUSIVE;
        // cout << "Allocate after BUS_RD at block 0x" << hex << trans.block_addr << dec
        //      << ": state set to " << (trans.found_in_other ? "SHARED" : "EXCLUSIVE") << endl;
      } else if (trans.op == BUS_RDX) {
        line.state = MODIFIED;
        line.dirty = true;
        // cout << "Allocate after BUS_RDX at block 0x" << hex << trans.block_addr << dec
        //      << ": state set to MODIFIED" << endl;
      }
      update_lru(set_index, way);
    }
  }

  // **Check if cache has a block**
  bool has_block(uint32_t block_addr) const {
    uint32_t set_index = (block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = block_addr >> (block_offset_bits + set_index_bits);
    for (const auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) return true;
    }
    return false;
  }
  MESIState get_line_state(uint32_t block_addr) const {
    // decode set index and tag
    uint32_t set_index = (block_addr >> block_offset_bits) & ((1u << set_index_bits) - 1);
    uint32_t tag = block_addr >> (block_offset_bits + set_index_bits);

    // scan all ways in the set
    for (const auto& line : cache_array[set_index]) {
      if (line.state != INVALID && line.tag == tag) {
        return line.state;
      }
    }

    // not found
    return INVALID;
  }

  int get_miss_count() const { return miss_count; }
  int get_access_count() const { return access_count; }
  int get_eviction_count() const { return eviction_count; }
  int get_writeback_count() const { return writeback_count; }
  int get_block_size() const { return block_size; }
};

// **Bus Class: Manages bus transactions with randomized order**
class Bus {
 private:
  BusTransaction current_transaction;  // Current active transaction
  int remaining_cycles;                // Cycles left for current transaction
  deque<BusTransaction> pending;       // queue of pending transactions
  vector<BusTransaction> cycle_transactions;
  int invalidations;      // Total number of invalidations
  uint64_t data_traffic;  // Total bus traffic in bytes
  mt19937 rng;            // Random number generator for shuffling

 public:
  unordered_set<uint32_t> pending_read_blocks;      // blocks with an outstanding BUS_RD
  unordered_map<uint32_t, vector<int>> rd_waiters;  // cores waiting on each block  vector<BusTransaction> cycle_transactions;  // Transactions collected in current cycle

  Bus() : remaining_cycles(0), invalidations(0), data_traffic(0), current_transaction() {
    random_device rd;
    rng = mt19937(rd());
  }

  // **Collect transaction for the current cycle**
  void add_transaction(const BusTransaction& t) {
    if (t.op == BUS_RD) {
      // If a read to this block is already pending, register the core as a waiter
      auto block = t.block_addr;
      if (pending_read_blocks.count(block)) {
        rd_waiters[block].push_back(t.core_id);
        return;  // skip enqueue
      }
      // First read for this block: mark it and enqueue
      pending_read_blocks.insert(block);
    }
    // Collect for this cycle
    cycle_transactions.push_back(t);
  }

  // **Randomize and queue transactions collected in the current cycle**
  void commit_transactions() {
    for (auto& t : cycle_transactions) {
      pending.push_back(t);  // now using deque
    }
    cycle_transactions.clear();
  }

  // **Start the next transaction**
  void start_next_transaction(const std::vector<Cache*>& caches) {
    if (pending.empty()) return;
    current_transaction = pending.front();
    pending.pop_front();

    uint32_t block = current_transaction.block_addr;
    bool found = false;
    bool foundM = false;
    // Scan all other caches for a copy, and note if any is MODIFIED
    for (int i = 0; i < (int)caches.size(); i++) {
      if (i == current_transaction.core_id) continue;
      if (caches[i]->has_block(block)) {
        found = true;
        // look up the line’s state
        MESIState st = caches[i]->get_line_state(block);
        if (st == MODIFIED) {
          foundM = true;
          break;  // a modified copy dominates
        }
      }
    }
    current_transaction.found_in_other = found;

    int block_size = caches[0]->get_block_size();
    int N = block_size / 4;  // e.g. 64-byte block → N=16

    if (current_transaction.op == BUS_RD) {
      if (!found) {
        // (1) No other copies → memory fetch
        remaining_cycles = 100;
      } else {
        // (2,3,4) Someone has it: abort memory, do cache-to-cache
        // If it was MODIFIED, they must write back first (100 cycles)
        remaining_cycles = (foundM ? 100 /* writeback of M */ : 0) + 2 * N;  // cache-to-cache transfer
      }
    } else {  // BUS_RDX (write-miss) – unchanged from before
      if (found) {
        remaining_cycles = 100 + (foundM ? 100 : 0);
      } else {
        remaining_cycles = 100;
      }
      if (current_transaction.needs_writeback)
        remaining_cycles += 100;
    }

    data_traffic += block_size;
  }

  // **Process one cycle of the transaction**
  bool process_cycle() {
    if (remaining_cycles > 0) {
      remaining_cycles--;
      if (remaining_cycles == 0) {
        // cout << "Bus transaction completed: core " << current_transaction.core_id
        //      << ", op " << (current_transaction.op == BUS_RD ? "BUS_RD" : "BUS_RDX")
        //      << ", address 0x" << hex << current_transaction.address << dec << endl;
      }
      return remaining_cycles == 0;
    }
    return false;
  }

  bool is_busy() const { return remaining_cycles > 0; }
  BusTransaction& get_current_transaction() { return current_transaction; }
  int get_invalidations() const { return invalidations; }
  uint64_t get_data_traffic() const { return data_traffic; }
  void increment_invalidations() { invalidations++; }
  bool is_pending_empty() const { return pending.empty(); }
  bool has_pending() const { return !pending.empty(); }
};

// **Core Class: Represents a processor core**
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

  bool has_next_instruction() { return trace_file.peek() != EOF; }

  // **Process one instruction**
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
  void increment_idle_cycles() { idle_cycles++, total_cycles++; }
  int get_id() const { return id; }
  int get_instruction_count() const { return instruction_count; }
  uint64_t get_total_cycles() const { return total_cycles; }
  void increment_total_cycles() { total_cycles++; }
  uint64_t get_idle_cycles() const { return idle_cycles; }
  double get_miss_rate() const {
    return cache->get_access_count() > 0 ? static_cast<double>(cache->get_miss_count()) / cache->get_access_count() : 0.0;
  }
  int get_eviction_count() const { return cache->get_eviction_count(); }
  int get_writeback_count() const { return cache->get_writeback_count(); }
  Cache* get_cache() { return cache; }
  bool get_is_stalled() const { return is_stalled; }
};

// **Simulator Class: Orchestrates the simulation**
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

  // **Run the simulation**
  void run() {
    while (true) {
      // Process bus transaction
      if (bus->is_busy()) {
        if (bus->process_cycle()) {
          auto& done = bus->get_current_transaction();
          auto block = done.block_addr;

          if (done.op == BUS_RD) {
            // Wake initiator
            cores[done.core_id]->set_stalled(false);
            cores[done.core_id]->get_cache()->update_after_bus(done);

            // Wake any additional waiters on this block
            auto it = bus->rd_waiters.find(block);
            if (it != bus->rd_waiters.end()) {
              for (int cid : it->second) {
                cores[cid]->set_stalled(false);
                cores[cid]->get_cache()->update_after_bus(done);
              }
              bus->rd_waiters.erase(it);
            }
          } else {
            // For BUS_RDX or other ops, handle as before
            cores[done.core_id]->set_stalled(false);
            cores[done.core_id]->get_cache()->update_after_bus(done);
          }

          // Standard snoop/invalidation for all caches
          for (int i = 0; i < 4; i++) {
            if (i != done.core_id && cores[i]->get_cache()->snoop(done)) {
              bus->increment_invalidations();
            }
          }

          // Remove the block from pending_read_blocks on completion
          if (done.op == BUS_RD) {
            bus->pending_read_blocks.erase(block);
          }
        }
      }

      // Start next transaction if bus is idle
      if (!bus->is_busy() && bus->has_pending()) {
        bus->start_next_transaction(get_caches());
      }

      // Process each core
      for (auto& core : cores) {
        int cid = core->get_id();

        // Case 1: core is actively processing
        if (!core->get_is_stalled() && core->has_next_instruction()) {
          core->increment_total_cycles();
          core->process_instruction(bus);
        }
        // Case 2: core is stalled
        else if (core->get_is_stalled()) {
          // Only increment idle_cycles if this core is *not* the one owning the bus
          if (!bus->is_busy() || cid != bus->get_current_transaction().core_id) {
            core->increment_idle_cycles();
          } else {
            core->increment_total_cycles();
          }
        }
        // Case 3: core finished its trace
        // Do nothing
      }

      // Commit transactions for this cycle
      bus->commit_transactions();

      // Check if simulation is complete
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

  // **Write output to file**
  void write_output() {
    ofstream out(output_file);
    for (int i = 0; i < 4; i++) {
      out << "Core " << i << ":\n";
      out << "Instructions: " << cores[i]->get_instruction_count() << "\n";
      out << "Total Cycles: " << cores[i]->get_total_cycles() << ":\n";
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

// **Command-line argument parsing**
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