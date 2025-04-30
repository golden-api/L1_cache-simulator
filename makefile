CC = g++
CFLAGS = -Wall -std=c++11 -O2
TARGET = L1simulate
SOURCES = q.cpp
OBJECTS = $(SOURCES:.cpp=.o)
LATEX = pdflatex
LATEXFLAGS = -interaction=nonstopmode -halt-on-error
REPORT = report.tex
PDF = report.pdf

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

report: $(PDF)

$(PDF): $(REPORT)
	@$(LATEX) $(LATEXFLAGS) $(REPORT) > /dev/null
	@$(LATEX) $(LATEXFLAGS) $(REPORT) > /dev/null

clean:
	rm -f $(OBJECTS) $(TARGET) $(PDF) *.aux *.log *.out

test: $(TARGET)
	./$(TARGET) -t app1 -s 6 -E 2 -b 5 -o output.txt

.PHONY: all clean test report
