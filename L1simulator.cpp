#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

// Enums for MESI states and operation types
enum MESIState { INVALID,
                 SHARED,
                 EXCLUSIVE,
                 MODIFIED };
enum OpType { READ,
              WRITE };
enum BusOp { BUS_RD,
             BUS_RDX,
             BUS_UPDATE };

// Structure to represent a cache line
struct CacheLine {
  uint32_t tag;
  MESIState state;
  bool dirty;
  int lru_order;  // Lower value means more recently used
  CacheLine() : tag(0), state(INVALID), dirty(false), lru_order(0) {}
};

// Structure for bus transactions
struct BusTransaction {
  int core_id;
  BusOp op;
  uint32_t address;
  uint32_t block_addr;  // Address of the block (tag + index)
  int cycles_remaining;
  BusTransaction(int id, BusOp o, uint32_t addr, uint32_t baddr, int cycles)
      : core_id(id), op(o), address(addr), block_addr(baddr), cycles_remaining(cycles) {}
};

// Class to manage a single L1 cache
class Cache {
 public:
  int num_sets;
  int assoc;
  int block_size;
  int set_index_bits;
  int block_offset_bits;
  std::vector<std::vector<CacheLine>> sets;
  int miss_count;
  int access_count;
  int eviction_count;
  int writeback_count;

 public:
  Cache(int s, int e, int b) : num_sets(1 << s), assoc(e), block_size(1 << b), set_index_bits(s), block_offset_bits(b), miss_count(0), access_count(0), eviction_count(0), writeback_count(0) {
    sets.resize(num_sets, std::vector<CacheLine>(assoc));
  }

  // Access the cache for a read or write
  bool access(uint32_t address, OpType op, int& cycles, bool& needs_bus, BusOp& bus_op, uint32_t& block_addr) {
    access_count++;
    uint32_t index = (address >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = address >> (set_index_bits + block_offset_bits);
    block_addr = address & ~(block_size - 1);  // Block-aligned address

    // Check for hit
    for (int i = 0; i < assoc; i++) {
      CacheLine& line = sets[index][i];
      if (line.state != INVALID && line.tag == tag) {
        // Hit
        update_lru(index, i);
        if (op == READ) {
          if (line.state == MODIFIED || line.state == EXCLUSIVE || line.state == SHARED) {
            cycles = 1;  // Cache hit
            return true;
          }
        } else {  // WRITE
          if (line.state == MODIFIED || line.state == EXCLUSIVE) {
            line.dirty = true;
            line.state = MODIFIED;
            cycles = 1;
            return true;
          } else if (line.state == SHARED) {
            // Need to invalidate other caches
            needs_bus = true;
            bus_op = BUS_RDX;
            line.dirty = true;
            line.state = MODIFIED;
            cycles = 1;
            return true;
          }
        }
      }
    }

    // Miss
    miss_count++;
    needs_bus = true;
    bus_op = (op == READ) ? BUS_RD : BUS_RDX;
    cycles = 100;  // Memory fetch
    int victim = find_lru_victim(index);
    CacheLine& line = sets[index][victim];

    // Handle eviction
    if (line.state != INVALID) {
      eviction_count++;
      if (line.state == MODIFIED) {
        writeback_count++;
        cycles += 100;  // Writeback to memory
      }
    }

    // Allocate new line
    line.tag = tag;
    line.state = (op == READ) ? EXCLUSIVE : MODIFIED;
    line.dirty = (op == WRITE);
    update_lru(index, victim);
    return false;
  }

  // Snoop bus transactions
  bool snoop(BusTransaction& trans, int& cycles, bool& provides_data, uint32_t& block_addr, int core_id) {
    uint32_t index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
    uint32_t tag = trans.block_addr >> (set_index_bits + block_offset_bits);
    provides_data = false;

    for (int i = 0; i < assoc; i++) {
      CacheLine& line = sets[index][i];
      if (line.state != INVALID && line.tag == tag) {
        if (trans.op == BUS_RD) {
          if (line.state == MODIFIED) {
            provides_data = true;
            cycles = block_size / 2;  // Transfer block
            line.state = SHARED;
            block_addr = trans.block_addr;
            return true;
          } else if (line.state == EXCLUSIVE) {
            line.state = SHARED;
            return true;
          } else if (line.state == SHARED) {
            return true;
          }
        } else if (trans.op == BUS_RDX) {
          if (line.state != INVALID) {
            line.state = INVALID;
            return true;
          }
        } else if (trans.op == BUS_UPDATE) {
          if (line.state == SHARED) {
            cycles = 2;  // Word transfer
            return true;
          }
        }
      }
    }
    return false;
  }

  // Update cache state after bus transaction
  void update_after_bus(BusTransaction& trans, bool got_data) {
    if (trans.op == BUS_RD || trans.op == BUS_RDX) {
      uint32_t index = (trans.block_addr >> block_offset_bits) & ((1 << set_index_bits) - 1);
      uint32_t tag = trans.block_addr >> (set_index_bits + block_offset_bits);
      for (int i = 0; i < assoc; i++) {
        CacheLine& line = sets[index][i];
        if (line.state != INVALID && line.tag == tag) {
          if (trans.op == BUS_RD && got_data) {
            line.state = SHARED;
          } else if (trans.op == BUS_RDX) {
            line.state = MODIFIED;
            line.dirty = true;
          }
          update_lru(index, i);
          break;
        }
      }
    }
  }

 private:
  // Update LRU order
  void update_lru(int index, int way) {
    for (int i = 0; i < assoc; i++) {
      if (sets[index][i].state != INVALID && i != way) {
        sets[index][i].lru_order++;
      }
    }
    sets[index][way].lru_order = 0;
  }

  // Find LRU victim
  int find_lru_victim(int index) {
    int max_lru = -1;
    int victim = 0;
    for (int i = 0; i < assoc; i++) {
      if (sets[index][i].state == INVALID) return i;
      if (sets[index][i].lru_order > max_lru) {
        max_lru = sets[index][i].lru_order;
        victim = i;
      }
    }
    return victim;
  }

 public:
  // Getters for statistics
  int get_miss_count() const { return miss_count; }
  int get_access_count() const { return access_count; }
  int get_eviction_count() const { return eviction_count; }
  int get_writeback_count() const { return writeback_count; }
  int get_block_size() const { return block_size; }
};

// Class to manage the bus
class Bus {
 private:
  std::queue<BusTransaction> pending;
  int invalidations;
  uint64_t data_traffic;

 public:
  Bus() : invalidations(0), data_traffic(0) {}

  void add_transaction(BusTransaction trans) {
    pending.push(trans);
  }

  bool process_transaction(const std::vector<Cache*>& caches, int& cycles, int& inv_count, uint64_t& traffic) {
    if (pending.empty()) return false;
    BusTransaction trans = pending.front();
    cycles = trans.cycles_remaining;
    bool got_data = false;
    inv_count = 0;

    for (int i = 0; i < caches.size(); i++) {
      if (i == trans.core_id) continue;
      int snoop_cycles = 0;
      bool provides_data = false;
      uint32_t block_addr = 0;
      bool hit = caches[i]->snoop(trans, snoop_cycles, provides_data, block_addr, i);
      if (hit) {
        if (trans.op == BUS_RDX && caches[i]->snoop(trans, snoop_cycles, provides_data, block_addr, i)) {
          inv_count++;
        }
        if (provides_data) {
          got_data = true;
          cycles = snoop_cycles;
          traffic += caches[i]->get_block_size();
        } else if (trans.op == BUS_UPDATE) {
          traffic += 4;  // Word transfer
        }
      }
    }

    if (!got_data && (trans.op == BUS_RD || trans.op == BUS_RDX)) {
      traffic += caches[trans.core_id]->get_block_size();  // From memory
    }

    if (trans.core_id != -1) {  // Ensure core_id is valid
      caches[trans.core_id]->update_after_bus(trans, got_data);
    }
    pending.pop();
    invalidations += inv_count;
    data_traffic += traffic;
    return true;
  }

  bool is_busy() const { return !pending.empty(); }

  int get_invalidations() const { return invalidations; }
  uint64_t get_data_traffic() const { return data_traffic; }
};

// Class to represent a core
class Core {
 private:
  std::ifstream trace_file;
  int instruction_count;
  uint64_t total_cycles;
  uint64_t idle_cycles;
  Cache* cache;
  bool blocked;
  int block_cycles_remaining;
  uint32_t current_addr;
  OpType current_op;

 public:
  Core(const std::string& trace_filename, int s, int e, int b) : instruction_count(0),
                                                                 total_cycles(0),
                                                                 idle_cycles(0),
                                                                 blocked(false),
                                                                 block_cycles_remaining(0) {
    cache = new Cache(s, e, b);
    trace_file.open(trace_filename);
    if (!trace_file.is_open()) {
      std::cerr << "Error opening trace file: " << trace_filename << std::endl;
      exit(1);
    }
  }

  ~Core() {
    if (trace_file.is_open()) trace_file.close();
    delete cache;
  }

  bool process_instruction(Bus* bus) {
    if (blocked) {
      block_cycles_remaining--;
      if (block_cycles_remaining <= 0) blocked = false;
      total_cycles++;
      idle_cycles++;
      return false;
    }

    if (!has_next_instruction()) return false;

    std::string line;
    std::getline(trace_file, line);
    std::istringstream iss(line);
    char op;
    uint32_t addr;
    iss >> op >> std::hex >> addr;
    current_op = (op == 'R') ? READ : WRITE;
    current_addr = addr;
    instruction_count++;

    int cycles = 0;
    bool needs_bus = false;
    BusOp bus_op;
    uint32_t block_addr;
    bool hit = cache->access(addr, current_op, cycles, needs_bus, bus_op, block_addr);

    total_cycles += cycles;
    if (cycles > 1) idle_cycles += (cycles - 1);

    if (needs_bus) {
      bus->add_transaction(BusTransaction(-1, bus_op, addr, block_addr, cycles));
      blocked = true;
      block_cycles_remaining = cycles;
    }

    return true;
  }

  void add_cycles(uint64_t cycles) {
    total_cycles += cycles;
    if (cycles > 0 && !blocked) idle_cycles += cycles;  // Adjust idle cycles as needed
  }

  bool has_next_instruction() {
    return trace_file.peek() != EOF;
  }

  // Getters
  int get_instruction_count() const { return instruction_count; }
  uint64_t get_total_cycles() const { return total_cycles; }
  uint64_t get_idle_cycles() const { return idle_cycles; }
  double get_miss_rate() const {
    return cache->get_access_count() > 0 ? static_cast<double>(cache->get_miss_count()) / cache->get_access_count() : 0.0;
  }
  int get_eviction_count() const { return cache->get_eviction_count(); }
  int get_writeback_count() const { return cache->get_writeback_count(); }
  Cache* get_cache() { return cache; }
};

// Class to manage the simulation
class Simulator {
 private:
  std::vector<Core*> cores;
  Bus* bus;
  std::string app_name;
  int set_index_bits;
  int assoc;
  int block_bits;
  std::string output_file;

 public:
  Simulator(const std::string& app, int s, int e, int b, const std::string& out)
      : app_name(app), set_index_bits(s), assoc(e), block_bits(b), output_file(out) {
    bus = new Bus();
    for (int i = 0; i < 4; i++) {
      std::string trace = app_name + "_proc" + std::to_string(i) + ".trace";
      cores.push_back(new Core(trace, s, e, b));
    }
  }

  ~Simulator() {
    for (Core* core : cores) delete core;
    delete bus;
  }

  void run() {
    bool active = true;
    while (active || bus->is_busy()) {
      active = false;
      int bus_cycles = 0;
      int inv_count = 0;
      uint64_t traffic = 0;

      // Process bus transactions first
      if (bus->process_transaction(get_caches(), bus_cycles, inv_count, traffic)) {
        for (Core* core : cores) {
          core->add_cycles(bus_cycles);
        }
      } else {
        // Process instructions
        for (Core* core : cores) {
          if (core->process_instruction(bus)) {
            active = true;
          } else {
            core->add_cycles(1);
          }
        }
      }
    }
  }

  void write_output() {
    std::ofstream out(output_file);
    if (!out.is_open()) {
      std::cerr << "Error opening output file: " << output_file << std::endl;
      return;
    }

    for (int i = 0; i < 4; i++) {
      out << "Core " << i << ":\n";
      out << "  Instructions: " << cores[i]->get_instruction_count() << "\n";
      out << "  Total Cycles: " << cores[i]->get_total_cycles() << "\n";
      out << "  Idle Cycles: " << cores[i]->get_idle_cycles() << "\n";
      out << "  Miss Rate: " << std::fixed << std::setprecision(4) << cores[i]->get_miss_rate() << "\n";
      out << "  Evictions: " << cores[i]->get_eviction_count() << "\n";
      out << "  Writebacks: " << cores[i]->get_writeback_count() << "\n";
    }
    out << "System:\n";
    out << "  Invalidations: " << bus->get_invalidations() << "\n";
    out << "  Bus Traffic (bytes): " << bus->get_data_traffic() << "\n";
    out.close();
  }

 private:
  std::vector<Cache*> get_caches() {
    std::vector<Cache*> caches;
    for (Core* core : cores) caches.push_back(core->get_cache());
    return caches;
  }
};

// Custom argument parser
struct Args {
  std::string app_name;
  int set_index_bits;
  int assoc;
  int block_bits;
  std::string output_file;
  bool show_help;
};

void print_help() {
  std::cout << "./L1simulate\n"
            << "-t <tracefile>: name of parallel application (e.g., app1)\n"
            << "-s <s>: number of set index bits\n"
            << "-E <E>: associativity\n"
            << "-b <b>: number of block bits\n"
            << "-o <outfilename>: output file\n"
            << "-h: prints this help\n";
}

Args parse_args(int argc, char* argv[]) {
  Args args = {"", 0, 0, 0, "", false};

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-h") {
      args.show_help = true;
      return args;
    }
    if (i + 1 >= argc) {
      std::cerr << "Missing value for option " << arg << std::endl;
      args.show_help = true;
      return args;
    }
    std::string value = argv[i + 1];

    if (arg == "-t") {
      args.app_name = value;
    } else if (arg == "-s") {
      try {
        args.set_index_bits = std::stoi(value);
      } catch (...) {
        std::cerr << "Invalid value for -s: " << value << std::endl;
        args.show_help = true;
      }
    } else if (arg == "-E") {
      try {
        args.assoc = std::stoi(value);
      } catch (...) {
        std::cerr << "Invalid value for -E: " << value << std::endl;
        args.show_help = true;
      }
    } else if (arg == "-b") {
      try {
        args.block_bits = std::stoi(value);
      } catch (...) {
        std::cerr << "Invalid value for -b: " << value << std::endl;
        args.show_help = true;
      }
    } else if (arg == "-o") {
      args.output_file = value;
    } else {
      std::cerr << "Unknown option: " << arg << std::endl;
      args.show_help = true;
    }
    ++i;  // Skip the value
  }

  return args;
}

int main(int argc, char* argv[]) {
  Args args = parse_args(argc, argv);

  if (args.show_help || args.app_name.empty() || args.set_index_bits <= 0 ||
      args.assoc <= 0 || args.block_bits <= 0 || args.output_file.empty()) {
    print_help();
    return 1;
  }

  Simulator sim(args.app_name, args.set_index_bits, args.assoc, args.block_bits, args.output_file);
  sim.run();
  sim.write_output();

  return 0;
}