# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Executable name
TARGET = L1simulate

# Source files
SRCS = main.cpp simulator.cpp core.cpp cache.cpp mesi.cpp bus.cpp
OBJS = $(SRCS:.cpp=.o)

# Header dependencies (if you include common headers)
DEPS = simulator.hpp core.hpp cache.hpp mesi.hpp bus.hpp

# Build the executable
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $<

# Clean up object files and executable
clean:
	rm -f $(OBJS) $(TARGET)

# Optional: run with example args
run:
	./L1simulate -t app1 -s 5 -E 4 -b 6 -o output.txt

.PHONY: all clean run
