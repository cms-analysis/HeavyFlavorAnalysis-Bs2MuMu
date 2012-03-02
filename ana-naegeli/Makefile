# Makefile for ana-producer

ifdef ROOTSYS
  ROOTCINT      = $(ROOTSYS)/bin/rootcint
  ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
  ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
  ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
else 
  ROOTCINT      = rootcint
  ROOTCFLAGS    = $(shell root-config --cflags)
  ROOTLIBS      = $(shell root-config --libs)
  ROOTGLIBS     = $(shell root-config --glibs)
endif

ifdef SCRAM_ARCH
  CXX         := $(shell scramv1 tool info cxxcompiler | grep CXX= | sed s/CXX=//)
else 
  ifdef CXXCOMPILER
    CXX         := $(CXXCOMPILER)
  else 
    CXX         := c++
  endif
endif

CXXFLAGS      = -g -O3 -Wall -fPIC -pipe $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)

OBJS = ana-produce.o

obj/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

all: bin obj bin/ana-produce
	@true

bin/ana-produce: $(addprefix obj/,$(OBJS))
	$(CXX) -o $@ $(addprefix obj/,$(OBJS)) $(LIBS)

bin:
	@mkdir -p bin

obj:
	@mkdir -p obj

clean:
	rm -rf obj

distclean: clean
	rm -rf bin
