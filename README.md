# Cache Simulator for Quad-Core Processors with MESI Coherence

This project simulates an L1 data cache for a quad-core processor system, implementing the MESI cache coherence protocol. The simulator is written in C++ and is designed to handle memory traces from parallel applications running on four cores. Use Git to manage the project and commit regularly. Submit the repository (without memory traces) on Moodle as `entrynum1_entrynum2_assign3.zip`.

## Project Structure

- **Source Code**: `q.cpp` - The main simulator code.
- **Makefile**: Automates compilation, testing, and PDF generation.
- **LaTeX Report**: `report.tex` - LaTeX source for the project report.
- **Traces**: (Not included) Memory trace files for each core, e.g., `app1_proc0.trace`, `app1_proc1.trace`, etc.
- **Output**: Simulation results are written to a specified output file and printed to the console.

## Compilation

To compile the simulator, use the provided `Makefile`. Ensure you have `g++` installed.

```bash
make
```

This will generate the executable `L1simulate`.

### Generating the Report PDF

To compile the LaTeX report (`report.tex`) into `report.pdf`, ensure `pdflatex` is installed and run:

```bash
make report
```

Alternatively, `make all` will compile both the simulator and the PDF.

### Cleaning Build Artifacts

To remove compiled objects, the executable, and LaTeX auxiliary files (including `report.pdf`):

```bash
make clean
```

## Running the Simulator

The simulator requires command-line arguments to specify the trace files, cache parameters, and output file.

### Command-Line Arguments

- `-t <tracefile>`: Prefix of the trace files (e.g., `app1` for `app1_proc0.trace`, etc.).
- `-s <s>`: Number of set index bits (number of sets = 2^s).
- `-E <E>`: Associativity (number of cache lines per set).
- `-b <b>`: Number of block bits (block size = 2^b bytes).
- `-o <outfilename>`: Output file for simulation results.
- `-h`: Print help message.

### Example Usage

```bash
./L1simulate -t app1 -s 6 -E 2 -b 5 -o output.txt
```

This command simulates a 4KB 2-way set-associative cache with 32-byte blocks using traces prefixed with `app1` and writes results to `output.txt`.

### Testing

A `test` target is provided in the `Makefile` for quick testing with default parameters:

```bash
make test
```

This runs the simulator with `s=6`, `E=2`, `b=5` (4KB cache, 2-way associative, 32-byte blocks) on `app1` traces.

## Input Format

- **Trace Files**: Each core's trace file (e.g., `app1_proc0.trace`) contains lines like:

  ```
  R 0x7e1afe78
  W 0x7e1ac04c
  ```

  where `R` denotes a read and `W` denotes a write, followed by a 32-bit hexadecimal address.

## Output

The simulator outputs the following metrics for each core:

- Total instructions
- Total reads
- Total writes
- Total cycles
- Idle cycles
- Misses
- Miss rate
- Evictions
- Writebacks
- Invalidations
- Data traffic (bytes)

Additionally, bus-related metrics are provided:

- Bus transactions
- Bus traffic (bytes)

Results are both printed to the console and saved to the specified output file.

## Assumptions

- **Memory Addresses**: 32-bit addresses. If fewer than 8 hex digits, assume leading zeros.
- **Word Size**: 4 bytes.
- **Cache Policy**: Write-back, write-allocate, LRU replacement.
- **Coherence Protocol**: MESI.
- **Bus Arbitration**: Cores with the earliest `currentTime` are prioritized.
- **Cache Blocking**: On a miss, the core halts until the miss is resolved.

## Bonus: Generating Interesting Traces

For bonus marks, small hand-generated traces can be created to demonstrate phenomena like false sharing. For example:

- **False Sharing Trace**:
  - Core 0 writes to address `0x0000` (block 0).
  - Core 1 writes to address `0x0010` (same block).
  - Repeated reads/writes to these addresses to observe many invalidations and bus traffic.

Example trace files:

- `app5_proc0.trace`:

  ```
   R 0x00000000  
   W 0x00000004  
   R 0x00000008  
   W 0x0000000C  
  ```

- `app5_proc1.trace`:

  ```
  R 0x00000010
  W 0x00000014
  R 0x00000000  
  W 0x00000018
  ```

Run with:

```bash
./L1simulate -t app5 -s 6 -E 2 -b 5 -o false_sharing_output.txt
```

## Notes

- Ensure trace files are in the same directory as the executable or provide the correct path.
- The simulator assumes trace files are correctly formatted; invalid entries are skipped.
- For varying cache parameters in experiments, modify the command-line arguments accordingly.
