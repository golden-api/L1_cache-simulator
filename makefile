CC = g++
CFLAGS = -Wall -std=c++11 -O2
TARGET = L1simulate
SOURCES = q.cpp
OBJECTS = $(SOURCES:.cpp=.o)
LATEX = pdflatex
LATEXFLAGS = -interaction=nonstopmode -halt-on-error
REPORT = report.tex
PDF = report.pdf
APPS = app1 app2 app3 app4 app5

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
	rm -f $(OBJECTS) $(TARGET) $(PDF) *.aux *.log *.out out*.txt

test: $(TARGET)
	@i=1; \
	for app in $(APPS); do \
		echo "Running $$app -> out$$i.txt"; \
		./$(TARGET) -t $$app -s 6 -E 2 -b 5 -o out$$i.txt || exit 1; \
		i=$$((i + 1)); \
	done

.PHONY: all clean test report
