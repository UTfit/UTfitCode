ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
CXX           = g++
CXXFLAGS      = -g -Wall -fPIC
CXXFLAGS      = -g -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared
##SOFLAGS       += -static

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
NGLIBS        += -lASImage
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
NGLIBS        += -lASImage
GLIBS          = $(filter-out -lNew, $(NGLIBS))

.SUFFIXES: .cc,.C

.C.o:
	$(CXX) $(CXXFLAGS) -c -I. $<



CKM_fit:  analysis.cc  
##test:  test.cc  
# -----------------------------------------------------------------------------
	$(CXX) $(CXXFLAGS) -c -I. $<
	$(LD) $(LDFLAGS) -o 2Dplot analysis.o  $(GLIBS)

 
# ================================================================================
clean:
	rm -f 2Dplot analysis.o
	rm -f *~
