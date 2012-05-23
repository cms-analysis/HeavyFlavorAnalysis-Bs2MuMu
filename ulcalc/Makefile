###############
# ENVIRONMENT #
###############

ifdef ROOTSYS
	ROOTCFLAGS	= $(shell $(ROOTSYS)/bin/root-config --cflags)
	ROOTLIBS	= $(shell $(ROOTSYS)/bin/root-config --libs)
else
	ROOTCFLAGS	= $(shell root-config --cflags)
	ROOTLIBS	= $(shell root-config --libs)
endif

CXXFLAGS	= -g -O3 -Wall -Werror -fPIC -pipe
CXX			= g++
LDFLAGS		= -g -lRooFitCore -lRooFit -lRooStats -lFoam -lMinuit -lMathMore
LD			= $(CXX)

#########
# FILES #
#########

HEADERS = external_constants.h ul_estimate.h bplus_estimator.h bmm_estimator.h
ULCALC_OBJS = ul_main.o ul_estimate.o external_constants.o
ULEXT_OBJS = ulc_ext.o bplus_estimator.o bmm_estimator.o external_constants.o
GRID_SEARCH_OBJS = grid_search.o external_constants.o ul_estimate.o bmm_estimator.o
GEN_OBJS = tree_gen.o


#################
# PATTERN RULES #
#################

obj/%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@


###########
# TARGETS #
###########

# Currently we only build the upper limit application
ulcalc: bin/ulcalc
	@true

all: ulcalc ulc_ext grid_search tree_gen
	@true

ulc_ext: bin/ulc_ext
	@true

grid_search: bin/grid_search
	@true

tree_gen: bin/tree_gen
	@true

lib/libNCRootUtils.so:
	ln -s /Users/cn/Documents/PSI/Projects/rootutils/lib/libNCRootUtils.so lib/libNCRootUtils.so

bin/ulcalc: $(addprefix obj/,$(ULCALC_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(ULCALC_OBJS)) $(ROOTLIBS) $(LDFLAGS)

bin/grid_search: $(addprefix obj/,$(GRID_SEARCH_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(GRID_SEARCH_OBJS)) $(ROOTLIBS) $(LDFLAGS)

bin/ulc_ext: $(addprefix obj/,$(ULEXT_OBJS)) lib/libNCRootUtils.so
	$(LD) -o $@ $(addprefix OBJ/,$(ULEXT_OBJS)) $(ROOTLIBS) $(LDFLAGS) lib/libNCRootUtils.so

bin/tree_gen: $(addprefix obj/,$(GEN_OBJS)) lib/libNCRootUtils.so
	$(LD) -o $@ $(addprefix obj/,$(GEN_OBJS)) $(ROOTLIBS) $(LDFLAGS) lib/libNCRootUtils.so

clean:
	rm -rf obj/*

cleanall: clean
	rm -rf bin/*

distclean: cleanall
	@true
