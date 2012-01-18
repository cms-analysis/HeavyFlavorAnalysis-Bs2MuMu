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
LDFLAGS		= -g -lRooFitCore -lRooFit -lRooStats -lFoam -Llib -lMinuit -lMathMore
LD			= $(CXX)

#########
# FILES #
#########

HEADERS = external_constants.h ul_estimate.h bplus_estimator.h bmm_estimator.h
ULCALC_OBJS = ul_main.o ul_estimate.o external_constants.o
BPLUS_OBJS = bplus_main.o bplus_estimator.o external_constants.o
BMM_OBJS = bmm_main.o bmm_estimator.o external_constants.o
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

all: ulcalc bplus_est bmm_est grid_search tree_gen
	@true

bplus_est: bin/bplus_est
	@true

bmm_est: bin/bmm_est
	@true

grid_search: bin/grid_search
	@true

tree_gen: bin/tree_gen
	@true

bin/ulcalc: $(addprefix obj/,$(ULCALC_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(ULCALC_OBJS)) $(ROOTLIBS) $(LDFLAGS)

bin/grid_search: $(addprefix obj/,$(GRID_SEARCH_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(GRID_SEARCH_OBJS)) $(ROOTLIBS) $(LDFLAGS)

bin/bplus_est: $(addprefix obj/,$(BPLUS_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(BPLUS_OBJS)) $(ROOTLIBS) $(LDFLAGS) ../rootutils/lib/libNCRootUtils.so

bin/bmm_est: $(addprefix obj/,$(BMM_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(BMM_OBJS)) $(ROOTLIBS) $(LDFLAGS)

bin/tree_gen: $(addprefix obj/,$(GEN_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(GEN_OBJS)) $(ROOTLIBS) $(LDFLAGS) ../rootutils/lib/libNCRootUtils.so

clean:
	rm -rf obj/*

cleanall: clean
	rm -rf bin/*

distclean: cleanall
	@true
