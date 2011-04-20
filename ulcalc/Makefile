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
LDFLAGS		= -g -lRooFitCore -lRooStats -lFoam -Llib -lMinuit
LD			= $(CXX)

#########
# FILES #
#########

ULCALC_HEADERS = external_constants.h ul_estimate.h
ULCALC_OBJS = ul_main.o ul_estimate.o external_constants.o
#GRID_SEARCH_HEADERS = external_constants.h ul_estimate.h bmm_estimator.h bplus_estimator.h
GRID_SEARCH_OBJS = grid_search.o external_constants.o ul_estimate.o bmm_estimator.o bplus_estimator.o


#################
# PATTERN RULES #
#################

obj/%.o : %.cpp $(ULCALC_HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@


###########
# TARGETS #
###########

# Currently we only build the grid search application.
all: bin/ulcalc
	@true

bin/ulcalc: $(addprefix obj/,$(ULCALC_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(ULCALC_OBJS)) $(ROOTLIBS) $(LDFLAGS)

bin/grid_search: $(addprefix obj/,$(GRID_SEARCH_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(GRID_SEARCH_OBJS)) $(ROOTLIBS) $(LDFLAGS)

clean:
	rm -rf obj/*

cleanall: clean
	rm -rf bin/*
