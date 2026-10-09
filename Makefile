#Makefile

CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
ROOTFLAGS = $(shell root-config --cflags --libs)

ARCH = -arch x86_64

CXXFLAGS += $(ROOTFLAGS)
#CXXFLAGS += $(ARCH)

Run: Run.C CorTreeAnalyzer.o EventClass.o
	$(CXX) $(CXXFLAGS) -o Run Run.C CorTreeAnalyzer.o EventClass.o libStRefMultCorr.so

CorTreeAnalyzer.o: CorTreeAnalyzer.cxx
	$(CXX) $(CXXFLAGS) -c CorTreeAnalyzer.cxx

EventClass.o: EventClass.cxx
	$(CXX) $(CXXFLAGS) -c EventClass.cxx

CorTreeAnalyzer.o: CorTreeAnalyzer.h

EventClass.o: EventClass.h

.PHONY: clean
clean:
	rm -f Run Run.o CorTreeAnalyzer.o EventClass.o
