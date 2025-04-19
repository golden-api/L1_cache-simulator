CC = g++
CFLAGS = -std=c++11 -Wall -O2
TARGET = L1simulator

all: $(TARGET)

$(TARGET): L1simulator.cpp
	$(CC) $(CFLAGS) -o $(TARGET) L1simulator.cpp

clean:
	rm -f $(TARGET)