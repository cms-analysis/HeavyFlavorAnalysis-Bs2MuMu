###############
# ENVIRONMENT #
###############

ifdef ROOTSYS
	ROOTCINT	= $(ROOTSYS)/bin/rootcint
	ROOTCFLAGS	= $(shell $(ROOTSYS)/bin/root-config --cflags)
	ROOTLIBS	= $(shell $(ROOTSYS)/bin/root-config --libs)
else
	ROOTCINT	= rootcint
	ROOTCFLAGS	= $(shell root-config --cflags)
	ROOTLIBS	= $(shell root-config --libs)
endif

CXXFLAGS	= -g -O3 -Wall -Werror -fPIC -pipe
CXX		= g++
LDFLAGS		= -g -lRooFitCore -lRooFit -lRooStats -lFoam -lMinuit -lMathMore
SOFLAGS		= -dynamiclib -fPIC
LD		= $(CXX)

#########
# FILES #
#########

HEADERS = external_constants.h ul_estimate.h bplus_estimator.h bmm_estimator.h
ULCALC_OBJS = ul_main.o ul_estimate.o external_constants.o
ULEXT_OBJS = ulc_ext.o bplus_estimator.o bmm_estimator.o external_constants.o
GEN_OBJS = tree_gen.o


#################
# PATTERN RULES #
#################

obj/%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

%Dict.cpp : %.h
	$(ROOTCINT) -f $@ -c $<

###########
# TARGETS #
###########

# Currently we only build the upper limit application
ulcalc: bin/ulcalc
	@true

all: ulcalc ulc_ext tree_gen
	@true

ulc_ext: bin/ulc_ext
	@true

tree_gen: bin/tree_gen
	@true

lib/libNCRootUtils.so:
	rm -f lib/libNCRootUtils.so
	cd lib && ln -s ../../rootutils/lib/libNCRootUtils.so && cd -

bin/ulcalc: obj bin $(addprefix obj/,$(ULCALC_OBJS))
	$(LD) -o $@ $(addprefix obj/,$(ULCALC_OBJS)) $(ROOTLIBS) $(LDFLAGS)

bin/ulc_ext: obj bin $(addprefix obj/,$(ULEXT_OBJS)) lib/libNCRootUtils.so
	$(LD) -o $@ $(addprefix obj/,$(ULEXT_OBJS)) $(ROOTLIBS) $(LDFLAGS) lib/libNCRootUtils.so

bin/tree_gen: obj bin $(addprefix obj/,$(GEN_OBJS)) lib/libNCRootUtils.so
	$(LD) -o $@ $(addprefix obj/,$(GEN_OBJS)) $(ROOTLIBS) $(LDFLAGS) lib/libNCRootUtils.so

obj:
	mkdir -p obj

bin:
	mkdir -p bin

clean:
	rm -rf obj
	rm -f *Dict*

cleanall: clean
	rm -rf bin

distclean: cleanall
	@true
