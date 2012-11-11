#############
# VARIABLES #
#############

CXX = g++
CXXFLAGS = -g -Wall -Os -I/opt/local/include -fpic
LDFLAGS = -shared -lTreePlayer

ifdef ROOTSYS
	ROOTCINT	= $(ROOTSYS)/bin/rootcint
	ROOTCFLAGS	= $(shell $(ROOTSYS)/bin/root-config --cflags)
	ROOTLIBS	= $(shell $(ROOTSYS)/bin/root-config --libs)
else
	ROOTCINT	= rootcint
	ROOTCFLAGS	= $(shell root-config --cflags)
	ROOTLIBS	= $(shell root-config --libs)
endif

#########
# FILES #
#########

UTILS = NCRootUtils.o NCFunctions.o NCMatrix.o
UTILS_HEADER = ${UTILS:.o=.h}
UTILS_DICTS = ${UTILS:.o=Dict.o}

#################
# PATTERN RULES #
#################

%Dict.cpp : %.h %LinkDef.h
	$(ROOTCINT) -f $@ -c $^

%Dict.cpp : %.h
	$(ROOTCINT) -f $@ -c $^

obj/%.o: %.cpp $(UTILS_HEADER) obj
	$(CXX) -c $(ROOTCFLAGS) $(CXXFLAGS) $< -o $@

###########
# TARGETS #
###########

all: lib/libNCRootUtils.so
	@true

lib/libNCRootUtils.so: $(addprefix obj/,$(UTILS) $(UTILS_DICTS)) lib
	$(CXX) $(LDFLAGS) $(ROOTLIBS) $(addprefix obj/,$(UTILS) $(UTILS_DICTS)) -o lib/libNCRootUtils.so

clean:
	rm -rf obj lib rm *Dict.*

obj:
	mkdir -p obj

lib:
	mkdir -p lib

# Dummy target for Xcode
dicts: $(patsubst %.o,%.cpp,$(UTILS_DICTS)) ;
