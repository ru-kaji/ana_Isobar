#Makefile

CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-unused-private-field
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

ARCH = -arch x86_64

CXXFLAGS += $(ROOTFLAGS)
#CXXFLAGS += $(ARCH)

Run: Run.C CorTreeAnalyzer.o EventClass.o
	$(CXX) $(CXXFLAGS) -o Run Run.C CorTreeAnalyzer.o EventClass.o libStRefMultCorr.so $(ROOTLIBS)

CorTreeAnalyzer.o: CorTreeAnalyzer.cxx
	$(CXX) $(CXXFLAGS) -c CorTreeAnalyzer.cxx

EventClass.o: EventClass.cxx
	$(CXX) $(CXXFLAGS) -c EventClass.cxx

CorTreeAnalyzer.o: CorTreeAnalyzer.h

EventClass.o: EventClass.h

.PHONY: clean
clean:
	rm -f Run Run.o CorTreeAnalyzer.o EventClass.o
