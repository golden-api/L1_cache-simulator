#include <bits/stdc++.h>
using namespace std;

enum State { I,
             S,
             E,
             M };  // MESI states

struct CacheLine {
  uint32_t tag;
  State state;
  long long lastUsed;
  CacheLine() : tag(0), state(I), lastUsed(0) {}
};

class Cache {
 public:
  int numSets, assoc, blockSize, setIndexBits;
  vector<vector<CacheLine>> sets;
  long long useCounter;
  Cache(int s, int E, int b)
      : numSets(1 << s), assoc(E), blockSize(1 << b), setIndexBits(s), useCounter(0) {
    sets.resize(numSets, vector<CacheLine>(assoc));
  }
  // Find a valid line for blockId in the set (extract set index and tag correctly)
  int findLine(uint32_t blockId) {
    uint32_t setIdx = blockId & ((1 << setIndexBits) - 1);
    uint32_t tag = blockId >> setIndexBits;
    for (int i = 0; i < assoc; i++) {
      if (sets[setIdx][i].state != I &&
          sets[setIdx][i].tag == tag) {
        return setIdx * assoc + i;
      }
    }
    return -1;
  }
  // Allocate a line for blockId (evict LRU if needed)
  int allocateLine(uint32_t blockId, bool &evictedDirty, bool &lineEvicted) {
    // Compute set index from lower s bits of blockId
    uint32_t setIdx = blockId & ((1 << setIndexBits) - 1);
    uint32_t tag = blockId >> setIndexBits;
    // Initialize flags
    lineEvicted = false;
    evictedDirty = false;
    // Look for an invalid line first (no eviction needed)
    for (int i = 0; i < assoc; i++) {
      if (sets[setIdx][i].state == I) {
        // Found empty line: allocate without eviction
        sets[setIdx][i].tag = tag;
        sets[setIdx][i].state = I;
        return setIdx * assoc + i;
      }
    }
    // Evict LRU in this set (true eviction)
    int lruIdx = 0;
    long long minUsed = sets[setIdx][0].lastUsed;
    for (int i = 1; i < assoc; i++) {
      if (sets[setIdx][i].lastUsed < minUsed) {
        minUsed = sets[setIdx][i].lastUsed;
        lruIdx = i;
      }
    }
    // Mark eviction happened
    lineEvicted = true;
    evictedDirty = (sets[setIdx][lruIdx].state == M);
    // Store new tag, state reset to Invalid until filled by bus transaction
    sets[setIdx][lruIdx].tag = tag;
    sets[setIdx][lruIdx].state = I;
    return setIdx * assoc + lruIdx;
  }
  void touchLine(int idx) {
    int setIdx = idx / assoc;
    int lineIdx = idx % assoc;
    sets[setIdx][lineIdx].lastUsed = ++useCounter;
  }
  CacheLine &getLine(int idx) {
    int setIdx = idx / assoc;
    int lineIdx = idx % assoc;
    return sets[setIdx][lineIdx];
  }
};

struct Stats {
  long long instructions = 0, reads = 0, writes = 0;
  long long cycles = 0, idleCycles = 0;
  long long misses = 0, evictions = 0, writebacks = 0;
  long long invalidations = 0, dataTraffic = 0;
};

class Core;
class Bus {
 public:
  int numCores = 4;
  int blockSize;
  long long busNextFree = 0;
  long long transactions = 0;
  long long traffic = 0;
  vector<Cache *> caches;
  Bus(int block) : blockSize(block) {}
  void addCache(Cache *c) { caches.push_back(c); }

  // Handle BusRd: load block (Shared or Exclusive)
  long long handleBusRd(int cid, int allocIdx, uint32_t blockId,
                        long long requestTime, Core *cores[]);

  // Handle BusRdX: read with intent to modify (invalidate others)
  long long handleBusRdX(int cid, int allocIdx, uint32_t blockId,
                         long long requestTime, Core *cores[]);

  // Handle BusUpgr: invalidate shared copies (S->I) for a write hit in S state
  long long handleBusUpgr(int cid, uint32_t blockId, long long requestTime,
                          Core *cores[]);
};

class Core {
 public:
  int id;
  Cache cache;
  vector<pair<int, uint32_t>> trace;
  int pc = 0;
  Stats stats;
  long long currentTime = 0;
  Bus *bus;
  Core(int id, int s, int E, int b) : id(id), cache(s, E, b) {}
  bool hasNext() const { return pc < (int)trace.size(); }
};

vector<Core *> allCores;

// BusRd implementation (for read miss)
long long Bus::handleBusRd(int cid, int allocIdx, uint32_t blockId,
                           long long requestTime, Core *cores[]) {
  int ownerM = -1;
  vector<int> sharers;
  // Find if any other cache has the block (M, S, or E)
  for (int i = 0; i < numCores; i++) {
    if (i == cid) continue;
    Cache *c = caches[i];
    int idx = c->findLine(blockId);
    if (idx != -1) {
      State st = c->getLine(idx).state;
      if (st == M) {
        ownerM = i;
        break;
      } else if (st == S || st == E) {
        sharers.push_back(i);
      }
    }
  }
  long long start = max(requestTime, busNextFree);
  long long endTime = start;
  transactions++;

  if (ownerM != -1) {
    // Another core has it in M: write back to memory then share
    Cache *ownerCache = caches[ownerM];
    int ownerIdx = ownerCache->findLine(blockId);
    endTime = start + 100;  // write-back latency
    traffic += blockSize;
    cores[ownerM]->stats.writebacks++;
    // State transition: M -> S for owner (just read, so it becomes shared)
    ownerCache->getLine(ownerIdx).state = S;
    endTime += 2 * (blockSize / 4);  // send block (2 cycles/word)
    traffic += blockSize;
    // Requestor gets S (shared)
    caches[cid]->getLine(allocIdx).state = S;
    // Any E->S in other sharers
    for (int other : sharers) {
      Cache *c = caches[other];
      int idx = c->findLine(blockId);
      if (idx != -1 && c->getLine(idx).state == E) {
        c->getLine(idx).state = S;
      }
    }
    cores[cid]->stats.dataTraffic += blockSize;
  } else if (!sharers.empty()) {
    // Some cores have S or E: just get from one of them
    endTime = start + 2 * (blockSize / 4);
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = S;
    for (int other : sharers) {
      Cache *c = caches[other];
      int idx = c->findLine(blockId);
      if (idx != -1 && c->getLine(idx).state == E) {
        c->getLine(idx).state = S;
      }
    }
    cores[cid]->stats.dataTraffic += blockSize;
  } else {
    // No one has it: fetch from memory
    endTime = start + 100;
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = E;  // exclusive copy
    cores[cid]->stats.dataTraffic += blockSize;
  }

  busNextFree = endTime;
  return endTime;
}
long long Bus::handleBusRdX(int cid, int allocIdx, uint32_t blockId,
                            long long requestTime, Core *cores[]) {
  // Gather which other caches hold this block
  vector<int> sharers;
  int ownerM = -1;
  for (int i = 0; i < numCores; i++) {
    if (i == cid) continue;
    Cache *c = caches[i];
    int idx = c->findLine(blockId);
    if (idx != -1) {
      State st = c->getLine(idx).state;
      if (st == M) {
        ownerM = i;
      } else if (st == S || st == E) {
        sharers.push_back(i);
      }
    }
  }

  // Serialize on the bus
  long long start = max(requestTime, busNextFree);
  long long endTime = start;
  transactions++;

  // 1) If someone has Modified, they must write it back (100 cycles),
  //    then we read from memory (another 100).
  if (ownerM != -1) {
    Cache *ownerCache = caches[ownerM];
    int ownerIdx = ownerCache->findLine(blockId);
    // write-back dirty line to memory
    endTime = start + 100;
    traffic += blockSize;
    cores[ownerM]->stats.writebacks++;
    ownerCache->getLine(ownerIdx).state = I;
    // Invalidate owner's copy
    cores[cid]->stats.invalidations++;

    // now fetch from memory
    endTime += 100;
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = M;
  }
  // 2) If others have Shared or Exclusive, just invalidate them and read from memory
  else if (!sharers.empty()) {
    // invalidate all sharers
    for (int other : sharers) {
      Cache *c = caches[other];
      int idx = c->findLine(blockId);
      if (idx != -1) {
        c->getLine(idx).state = I;
        cores[cid]->stats.invalidations++;
      }
    }
    // then fetch from memory
    endTime = start + 100;
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = M;
  }
  // 3) No other copies: simple memory fetch
  else {
    endTime = start + 100;
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = M;
  }

  busNextFree = endTime;
  return endTime;
}

// BusUpgr implementation (for write hit in S: invalidate others)
long long Bus::handleBusUpgr(int cid, uint32_t blockId,
                             long long requestTime, Core *cores[]) {
  long long start = max(requestTime, busNextFree);
  long long endTime = start;
  transactions++;
  for (int i = 0; i < numCores; i++) {
    if (i == cid) continue;
    Cache *c = caches[i];
    int idx = c->findLine(blockId);
    if (idx != -1 && (c->getLine(idx).state == S || c->getLine(idx).state == E)) {
      c->getLine(idx).state = I;
      cores[cid]->stats.invalidations++;
    }
  }
  busNextFree = endTime;
  return endTime;
}

void printHelp() {
  cout << "Usage: ./L1simulate -t <tracefile> -s <s> -E <E> -b <b> -o <outfilename>\n";
  cout << "-t <tracefile>: name of parallel application (e.g., app1)\n";
  cout << "-s <s>: number of set index bits (number of sets = 2^s)\n";
  cout << "-E <E>: associativity (number of cache lines per set)\n";
  cout << "-b <b>: number of block bits (block size = 2^b)\n";
  cout << "-o <outfilename>: logs output in file\n";
  cout << "-h: prints this help\n";
}

int main(int argc, char *argv[]) {
  string tracePrefix, outFilename;
  int s = -1, E = -1, b = -1;
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-h") {
      printHelp();
      return 0;
    } else if (arg == "-t" && i + 1 < argc) {
      tracePrefix = argv[++i];
    } else if (arg == "-s" && i + 1 < argc) {
      s = stoi(argv[++i]);
    } else if (arg == "-E" && i + 1 < argc) {
      E = stoi(argv[++i]);
    } else if (arg == "-b" && i + 1 < argc) {
      b = stoi(argv[++i]);
    } else if (arg == "-o" && i + 1 < argc) {
      outFilename = argv[++i];
    }
  }
  if (tracePrefix.empty() || s < 0 || E < 0 || b < 0 || outFilename.empty()) {
    cerr << "Missing required arguments\n";
    printHelp();
    return 1;
  }
  int blockSize = 1 << b;

  Bus bus(blockSize);
  vector<Core *> cores;
  for (int i = 0; i < 4; i++) {
    Core *c = new Core(i, s, E, b);
    c->bus = &bus;
    cores.push_back(c);
    bus.addCache(&c->cache);
  }
  allCores = cores;

  // Load traces: e.g., app1_proc0.trace, app1_proc1.trace, ...
  for (int i = 0; i < 4; i++) {
    string filename = tracePrefix + "_proc" + to_string(i) + ".trace";
    ifstream in(filename);
    if (!in) {
      cerr << "Error opening " << filename << "\n";
      return 1;
    }
    string op;
    uint32_t addr;
    while (in >> op >> hex >> addr >> dec) {
      int label = (op == "R" ? 2 : op == "W" ? 3
                                             : -1);
      if (label == -1) continue;
      cores[i]->trace.emplace_back(label, addr);
    }
    in.close();
  }

  // Simulation
  while (true) {
    int nextCore = -1;
    long long nextTime = LLONG_MAX;
    for (int i = 0; i < 4; i++) {
      Core *c = cores[i];
      if (c->hasNext() && c->currentTime < nextTime) {
        nextTime = c->currentTime;
        nextCore = i;
      }
    }
    if (nextCore < 0) break;
    Core *c = cores[nextCore];

    int label = c->trace[c->pc].first;
    uint32_t addr = c->trace[c->pc].second;
    c->pc++;
    c->stats.instructions++;

    // No-op instruction (if any) takes 1 cycle
    if (label == 0) {
      c->currentTime += 1;
      c->stats.cycles += 1;
      continue;
    }

    uint32_t blockId = addr >> b;  // Compute block address (discard offset)
    if (label == 2) {              // Read
      c->stats.reads++;
      int idx = c->cache.findLine(blockId);
      if (idx != -1) {
        // Hit: 1 cycle
        c->stats.cycles += 1;
        c->currentTime += 1;
        c->cache.touchLine(idx);
      } else {
        // Miss
        c->stats.misses++;
        bool evDirty = false;
        bool evLine = false;
        int allocIdx = c->cache.allocateLine(blockId, evDirty, evLine);
        if (evLine) c->stats.evictions++;  // FIX: Count true eviction (replaced a valid block)
        // Handle dirty eviction writeback if needed
        if (evDirty) {
          long long flushStart = max(c->currentTime, bus.busNextFree);
          long long flushWait = flushStart - c->currentTime;
          if (flushWait > 0) {
            // Bus was busy: count these as idle and cycles
            c->stats.idleCycles += flushWait;
            c->stats.cycles += flushWait;  // include as core time (exception)
          }
          long long flushEnd = flushStart + 100;
          c->stats.cycles += 100;  // writing back takes 100 cycles
          bus.transactions++;
          bus.traffic += blockSize;
          c->stats.writebacks++;
          bus.busNextFree = flushEnd;
          c->currentTime = flushEnd;
        }
        // Issue BusRd
        long long startTime = c->currentTime;
        long long endTime = bus.handleBusRd(c->id, allocIdx, blockId,
                                            startTime, allCores.data());
        // Core was stalled waiting: count idle cycles (minus 1 for this instruction)
        c->stats.idleCycles += (endTime - startTime);
        // The memory operation itself does not count as core work beyond 1 cycle
        c->stats.cycles += 1;  // finishing the read
        c->currentTime = endTime;
        c->cache.touchLine(allocIdx);
      }
    } else if (label == 3) {  // Write
      c->stats.writes++;
      int idx = c->cache.findLine(blockId);
      if (idx != -1) {
        State st = c->cache.getLine(idx).state;
        if (st == M) {
          // Hit in M: just 1 cycle
          c->stats.cycles += 1;
          c->currentTime += 1;
          c->cache.touchLine(idx);
        } else if (st == E) {
          // Hit in E: upgrade to M, 1 cycle
          c->stats.cycles += 1;
          c->currentTime += 1;
          c->cache.getLine(idx).state = M;
          c->cache.touchLine(idx);
        } else if (st == S) {
          // Hit in S: broadcast invalidate (BusUpgr) then upgrade to M
          long long startTime = c->currentTime;
          long long endTime = bus.handleBusUpgr(c->id, blockId, startTime,
                                                allCores.data());
          if (endTime > startTime) {
            // Bus was busy with others
            c->stats.idleCycles += (endTime - startTime);
          }
          // After invalidation, perform write (1 cycle)
          c->stats.cycles += 1;
          c->currentTime = endTime + 1;
          c->cache.getLine(idx).state = M;
          c->cache.touchLine(idx);
        }
      } else {
        // Write miss
        c->stats.misses++;
        bool evDirty = false;
        bool evLine = false;
        int allocIdx = c->cache.allocateLine(blockId, evDirty, evLine);
        if (evLine) c->stats.evictions++;  // FIX: Count true eviction (replaced a valid block)
        // Handle dirty eviction writeback if needed
        if (evDirty) {
          long long flushStart = max(c->currentTime, bus.busNextFree);
          long long flushWait = flushStart - c->currentTime;
          if (flushWait > 0) {
            c->stats.idleCycles += flushWait;
            c->stats.cycles += flushWait;
          }
          long long flushEnd = flushStart + 100;
          c->stats.cycles += 100;
          bus.transactions++;
          bus.traffic += blockSize;
          c->stats.writebacks++;
          bus.busNextFree = flushEnd;
          c->currentTime = flushEnd;
        }
        // Issue BusRdX (RWITM): bring block and invalidate others
        long long startTime = c->currentTime;
        long long endTime = bus.handleBusRdX(c->id, allocIdx, blockId,
                                             startTime, allCores.data());
        c->stats.idleCycles += (endTime - startTime);
        c->stats.cycles += 1;  // complete the write
        c->currentTime = endTime;
        c->cache.getLine(allocIdx).state = M;
        c->cache.touchLine(allocIdx);
      }
    }
  }
  // Output to file and console
  ofstream out(outFilename);
  if (!out) {
    cerr << "Error opening output file " << outFilename << "\n";
    return 1;
  }
  for (int i = 0; i < 4; i++) {
    Core *c = cores[i];
    out << "Core " << i << ":\n";
    cout << "Core " << i << ":\n";
    out << "  Total instructions: " << c->stats.instructions << "\n";
    cout << "  Total instructions: " << c->stats.instructions << "\n";
    out << "  Total reads: " << c->stats.reads << "\n";
    cout << "  Total reads: " << c->stats.reads << "\n";
    out << "  Total writes: " << c->stats.writes << "\n";
    cout << "  Total writes: " << c->stats.writes << "\n";
    out << "  Total cycles: " << c->stats.cycles << "\n";
    cout << "  Total cycles: " << c->stats.cycles << "\n";
    out << "  Idle cycles: " << c->stats.idleCycles << "\n";
    cout << "  Idle cycles: " << c->stats.idleCycles << "\n";
    out << "  Misses: " << c->stats.misses << "\n";
    cout << "  Misses: " << c->stats.misses << "\n";
    double missRate = 0.0;
    if (c->stats.reads + c->stats.writes > 0) {
      missRate = (double)c->stats.misses / (c->stats.reads + c->stats.writes) * 100.0;
    }
    out << "  Miss rate: " << fixed << setprecision(2) << missRate << "%\n";
    cout << "  Miss rate: " << fixed << setprecision(2) << missRate << "%\n";
    out << "  Evictions: " << c->stats.evictions << "\n";
    cout << "  Evictions: " << c->stats.evictions << "\n";
    out << "  Writebacks: " << c->stats.writebacks << "\n";
    cout << "  Writebacks: " << c->stats.writebacks << "\n";
    out << "  Invalidations: " << c->stats.invalidations << "\n";
    cout << "  Invalidations: " << c->stats.invalidations << "\n";
    out << "  Data traffic: " << c->stats.dataTraffic << " bytes\n";
    cout << "  Data traffic: " << c->stats.dataTraffic << " bytes\n";
  }
  out << "Bus transactions: " << bus.transactions << "\n";
  cout << "Bus transactions: " << bus.transactions << "\n";
  out << "Bus traffic: " << bus.traffic << " bytes\n";
  cout << "Bus traffic: " << bus.traffic << " bytes\n";
  out.close();

  // Cleanup
  for (Core *c : cores) delete c;
  return 0;
}