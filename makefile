CC = g++
CFLAGS = -Wall -std=c++11 -O2
TARGET = L1simulate
SOURCES = q.cpp
OBJECTS = $(SOURCES:.cpp=.o)
LATEX = pdflatex
LATEXFLAGS = -interaction=nonstopmode
REPORT = report.tex
PDF = report.pdf

all: $(TARGET) $(PDF)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(PDF): $(REPORT)
	$(LATEX) $(LATEXFLAGS) $(REPORT)
	$(LATEX) $(LATEXFLAGS) $(REPORT)  # Run twice for references

clean:
	rm -f $(OBJECTS) $(TARGET) $(PDF) *.aux *.log *.out

test: $(TARGET)
	./$(TARGET) -t app1 -s 6 -E 2 -b 5 -o output.txt

.PHONY: all clean test