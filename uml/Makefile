ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -ansi -D_GNU_SOURCE -O2
LDFLAGS       = -g -lRooFit -lRooFitCore -lMinuit -lFoam -lGraf

SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

CXXFLAGS      += $(ROOTCFLAGS)

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR)

ifdef ROOFITSYS
CXX	         += -I$(ROOFITSYS)include
GLIBS        += -L$(ROOFITSYS)lib
endif

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ../../../../../lib/

###########
# TARGETS #
###########

main_pdf_choise: $(INCLUDEDIR)main_pdf_choise.cc pdf_analysis.o
		$(CXX) $(CXXFLAGS) -o main_pdf_choise.o pdf_analysis.o $(GLIBS)  $(LDFLAGS) $ $<

pdf_analysis.o: $(INCLUDEDIR)pdf_analysis.cpp
	$(CXX) $(CXXFLAGS) -c pdf_analysis.cpp  $ $<



clean:
	rm -f *.o *~


