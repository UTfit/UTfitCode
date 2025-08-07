ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
CXX           = g++
CXXFLAGS      = -g -Wall -fPIC
CXXFLAGS      = -g -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit 
NGLIBS        += -lRooFitCore
NGLIBS        += -lFoam
GLIBS          = $(filter-out -lNew, $(NGLIBS))

.SUFFIXES: .cc,.C

.C.o:
	$(CXX) $(CXXFLAGS) -c -I. $<



CKM_gen:  CKM_gen.cc  
##test:  test.cc
# -----------------------------------------------------------------------------
	$(CXX) $(CXXFLAGS) -c -I. $<
	$(LD) $(LDFLAGS) -o CKM_gen CKM_gen.o  $(GLIBS)

 
# ================================================================================
clean:
	rm -f CKM_gen.o  CKM_gen
	rm -f *~
