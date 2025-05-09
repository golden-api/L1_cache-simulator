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
  int findLine(uint32_t blockId) {
    uint32_t setIdx = blockId & ((1 << setIndexBits) - 1);
    uint32_t tag = blockId >> setIndexBits;
    for (int i = 0; i < assoc; i++) {
      if (sets[setIdx][i].state != I && sets[setIdx][i].tag == tag) {
        return setIdx * assoc + i;
      }
    }
    return -1;
  }
  int allocateLine(uint32_t blockId, bool &evictedDirty, bool &lineEvicted) {
    uint32_t setIdx = blockId & ((1 << setIndexBits) - 1);
    uint32_t tag = blockId >> setIndexBits;
    lineEvicted = false;
    evictedDirty = false;
    for (int i = 0; i < assoc; i++) {
      if (sets[setIdx][i].state == I) {
        sets[setIdx][i].tag = tag;
        sets[setIdx][i].state = I;
        return setIdx * assoc + i;
      }
    }
    int lruIdx = 0;
    long long minUsed = sets[setIdx][0].lastUsed;
    for (int i = 1; i < assoc; i++) {
      if (sets[setIdx][i].lastUsed < minUsed) {
        minUsed = sets[setIdx][i].lastUsed;
        lruIdx = i;
      }
    }
    lineEvicted = true;
    evictedDirty = (sets[setIdx][lruIdx].state == M);
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

  long long handleBusRd(int cid, int allocIdx, uint32_t blockId,
                        long long requestTime, Core *cores[]);
  long long handleBusRdX(int cid, int allocIdx, uint32_t blockId,
                         long long requestTime, Core *cores[]);
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
  for (int i = 0; i < numCores; i++) {
    if (i == cid) continue;
    Cache *c = caches[i];
    int idx = c->findLine(blockId);
    if (idx != -1) {
      State st = c->getLine(idx).state;
      if (st == M) {
        ownerM = i;
        break;
      } else if (st == S || st == E)
        sharers.push_back(i);
    }
  }
  long long start = max(requestTime, busNextFree);
  long long endTime = start;
  transactions++;

  if (ownerM != -1) {
    Cache *ownerCache = caches[ownerM];
    int ownerIdx = ownerCache->findLine(blockId);
    endTime = start + 100;
    traffic += blockSize;  // memory write-back
    cores[ownerM]->stats.writebacks++;
    ownerCache->getLine(ownerIdx).state = S;
    endTime += 2 * (blockSize / 4);
    traffic += blockSize;  // cache-to-cache transfer
    caches[cid]->getLine(allocIdx).state = S;
    for (int other : sharers) {
      Cache *c = caches[other];
      int idx = c->findLine(blockId);
      if (idx != -1 && c->getLine(idx).state == E) c->getLine(idx).state = S;
    }
    cores[cid]->stats.dataTraffic += blockSize;
  } else if (!sharers.empty()) {
    endTime = start + 2 * (blockSize / 4);
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = S;
    for (int other : sharers) {
      Cache *c = caches[other];
      int idx = c->findLine(blockId);
      if (idx != -1 && c->getLine(idx).state == E) c->getLine(idx).state = S;
    }
    cores[cid]->stats.dataTraffic += blockSize;
  } else {
    endTime = start + 100;
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = E;
    cores[cid]->stats.dataTraffic += blockSize;
  }

  busNextFree = endTime;
  return endTime;
}

long long Bus::handleBusRdX(int cid, int allocIdx, uint32_t blockId,
                            long long requestTime, Core *cores[]) {
  vector<int> sharers;
  int ownerM = -1;
  for (int i = 0; i < numCores; i++) {
    if (i == cid) continue;
    Cache *c = caches[i];
    int idx = c->findLine(blockId);
    if (idx != -1) {
      State st = c->getLine(idx).state;
      if (st == M)
        ownerM = i;
      else if (st == S || st == E)
        sharers.push_back(i);
    }
  }
  long long start = max(requestTime, busNextFree);
  long long endTime = start;
  transactions++;

  if (ownerM != -1) {
    Cache *ownerCache = caches[ownerM];
    int ownerIdx = ownerCache->findLine(blockId);
    endTime = start + 100;
    traffic += blockSize;  // write-back
    cores[ownerM]->stats.writebacks++;
    ownerCache->getLine(ownerIdx).state = I;
    cores[cid]->stats.invalidations++;
    endTime += 100;
    traffic += blockSize;  // memory fetch
    caches[cid]->getLine(allocIdx).state = M;
  } else if (!sharers.empty()) {
    for (int other : sharers) {
      Cache *c = caches[other];
      int idx = c->findLine(blockId);
      if (idx != -1) {
        c->getLine(idx).state = I;
        cores[cid]->stats.invalidations++;
      }
    }
    endTime = start + 100;
    traffic += blockSize;
    caches[cid]->getLine(allocIdx).state = M;
  } else {
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
  int num_sets = 1 << s;
  int cache_size = num_sets * E * blockSize;
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
    string filename = tracePrefix + "/" + tracePrefix + "_proc" + to_string(i) + ".trace";
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
  int lastCore = -1;
  // Simulation
  while (true) {
    int nextCore = -1;
    const int N = cores.size();
    for (int step = 1; step <= N; ++step) {
      int idx = (lastCore + step) % N;
      if (cores[idx]->hasNext()) {
        nextCore = idx;
        break;
      }
    }
    // If none have work, we’re done
    if (nextCore < 0) break;

    // 2. Update lastCore and dispatch
    lastCore = nextCore;
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
        c->stats.cycles += 1;
        c->stats.misses++;
        bool evDirty = false;
        bool evLine = false;
        int allocIdx = c->cache.allocateLine(blockId, evDirty, evLine);
        if (evLine) c->stats.evictions++;
        // Handle dirty eviction writeback if needed
        if (evDirty) {
          long long flushStart = max(c->currentTime, bus.busNextFree);
          long long flushWait = flushStart - c->currentTime;
          c->stats.idleCycles += flushWait;
          long long flushEnd = flushStart + 100;
          // removed bus.transactions++ and bus.traffic += blockSize;
          c->stats.writebacks++;
          bus.busNextFree = flushEnd;
          c->currentTime = flushEnd;
        }
        // Issue BusRd
        long long startTime = c->currentTime;
        long long endTime = bus.handleBusRd(c->id, allocIdx, blockId, startTime, allCores.data());
        c->stats.idleCycles += (endTime - startTime);
        c->stats.cycles += 1;
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
          long long endTime = bus.handleBusUpgr(c->id, blockId, startTime, allCores.data());
          if (endTime > startTime) {
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
        c->stats.cycles += 1;
        c->stats.misses++;
        bool evDirty = false;
        bool evLine = false;
        int allocIdx = c->cache.allocateLine(blockId, evDirty, evLine);
        if (evLine) c->stats.evictions++;
        if (evDirty) {
          long long flushStart = max(c->currentTime, bus.busNextFree);
          long long flushWait = flushStart - c->currentTime;
          c->stats.idleCycles += flushWait;
          long long flushEnd = flushStart + 100;
          c->stats.writebacks++;
          bus.busNextFree = flushEnd;
          c->currentTime = flushEnd;
        }
        // Issue BusRdX (RWITM): bring block and invalidate others
        long long startTime = c->currentTime;
        long long endTime = bus.handleBusRdX(c->id, allocIdx, blockId, startTime, allCores.data());
        c->stats.idleCycles += (endTime - startTime);
        c->stats.cycles += 1;
        c->currentTime = endTime;
        c->cache.getLine(allocIdx).state = M;
        c->cache.touchLine(allocIdx);
      }
    }
  }
  long long maxExecTime = 0;
  for (int i = 0; i < 4; ++i) {
    maxExecTime = max(maxExecTime, cores[i]->currentTime);
  }
  // Output to file and console
  ofstream out(outFilename);
  if (!out) {
    cerr << "Error opening output file " << outFilename << "\n";
    return 1;
  }
  out << "Simulation Parameters" << ":\n";
  out << "Trace Prefix: " << tracePrefix << "\n";
  out << "Set Index Bits: " << s << "\n";
  out << "Associativity: " << E << "\n";
  out << "Block Bits: " << b << "\n";
  out << "Block Size (Bytes): " << blockSize << "\n";
  out << "Number of Sets: " << num_sets << "\n";
  out << "Cache Size (KB per core): " << cache_size << "\n";
  out << "MESI Protocol: " << "Enabled" << "\n";
  out << "Write Policy: " << "Write-back" << "\n";
  out << "Replacement Policy: " << "LRU" << "\n";
  out << "Bus: " << "Central snooping bus" << "\n";
  out << "\n";

  for (int i = 0; i < 4; i++) {
    Core *c = cores[i];
    out << "Core " << i << " Statistics:\n";
    out << "Total instructions: " << c->stats.instructions << "\n";
    out << "Total reads: " << c->stats.reads << "\n";
    out << "Total writes: " << c->stats.writes << "\n";
    out << "Total Execution Cycles: " << c->stats.cycles << "\n";
    out << "Idle Cycles: " << c->stats.idleCycles << "\n";
    out << "Cache Misses: " << c->stats.misses << "\n";
    double missRate = 0.0;
    if (c->stats.reads + c->stats.writes > 0) {
      missRate = (double)c->stats.misses / (c->stats.reads + c->stats.writes) * 100.0;
    }
    out << "Cache Miss Rate: " << fixed << setprecision(2) << missRate << "%\n";
    out << "Cache Evictions: " << c->stats.evictions << "\n";
    out << "Writebacks: " << c->stats.writebacks << "\n";
    out << "Bus Invalidations: " << c->stats.invalidations << "\n";
    out << "Data traffic: " << c->stats.dataTraffic << " bytes\n";
    out << "\n";
    out << "Overall Bus Summary:\n";
  }
  out << "Total Bus Transactions: " << bus.transactions << endl;
  out << "Total Bus Traffic (Bytes): " << bus.traffic << endl;
  out.close();

  // Cleanup
  for (Core *c : cores) delete c;
  return 0;
}