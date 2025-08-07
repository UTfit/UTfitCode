ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
 
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
GLIBS          = $(filter-out -lNew, $(NGLIBS))

.SUFFIXES: .cc,.C

.C.o:
	$(CXX) $(CXXFLAGS) -c -I. $<



CKM_gen_rhopi:  CKM_gen_rhopi.cc  
##test:  test.cc
# -----------------------------------------------------------------------------
	$(CXX) $(CXXFLAGS) -c -I. $<
	$(LD) $(LDFLAGS) -o CKM_gen_rhopi CKM_gen_rhopi.o  $(GLIBS)

 
# ================================================================================
clean:
	rm -f CKM_gen_rhopi.o  CKM_gen_rhopi
	rm -f *~
