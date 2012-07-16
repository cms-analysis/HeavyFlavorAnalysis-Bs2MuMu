ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -fPIC -ansi -D_GNU_SOURCE -O2 -Wall
LDFLAGS       = -g -lRooFit -lRooFitCore -lMinuit -lFoam -lGraf -lRooStats

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

KSTest: $(INCLUDEDIR)KSTest.cc
	$(CXX) $(CXXFLAGS) -o KSTest.o $(GLIBS)  $(LDFLAGS) $ $<

main_toyMC: $(INCLUDEDIR)main_toyMC.cc pdf_analysis.o pdf_fitData.o pdf_toyMC.o
		$(CXX) $(CXXFLAGS) -o main_toyMC.o pdf_analysis.o pdf_fitData.o pdf_toyMC.o $(GLIBS)  $(LDFLAGS) $ $<

main_fitData: $(INCLUDEDIR)main_fitData.cc pdf_analysis.o pdf_fitData.o
		$(CXX) $(CXXFLAGS) -o main_fitData.o pdf_analysis.o pdf_fitData.o $(GLIBS)  $(LDFLAGS) $ $<

main_pdf_choise: $(INCLUDEDIR)main_pdf_choise.cc pdf_analysis.o
		$(CXX) $(CXXFLAGS) -o main_pdf_choise.o pdf_analysis.o $(GLIBS)  $(LDFLAGS) $ $<

pdf_analysis.o: $(INCLUDEDIR)pdf_analysis.cpp
	$(CXX) $(CXXFLAGS) -c pdf_analysis.cpp  $ $<

pdf_toyMC.o: $(INCLUDEDIR)pdf_toyMC.cpp
	$(CXX) $(CXXFLAGS) -c pdf_toyMC.cpp  $ $<

pdf_fitData.o: $(INCLUDEDIR)pdf_fitData.cpp
	$(CXX) $(CXXFLAGS) -c pdf_fitData.cpp  $ $<

clean:
	rm -f *.o *~


