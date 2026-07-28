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


# Two generators live side by side:
#
#   make            -> CKM_gen_rhopi       flat prior + importance weighting
#                                          (the original code, still the default)
#   make mcmc       -> CKM_gen_rhopi_mcmc  parallel tempering (replica-exchange
#                                          MCMC); needs the sampler block in the
#                                          .conf (NReplicas, NTemps, BetaMin,
#                                          NBurnin, SwapEvery, StepMag, StepPhase)
#   make all        -> both
#
# Both write the same 12 histograms into Input/, so downstream macros do not
# care which produced a given .root file.

.PHONY: all mcmc clean

.DEFAULT_GOAL := CKM_gen_rhopi

all: CKM_gen_rhopi CKM_gen_rhopi_mcmc

mcmc: CKM_gen_rhopi_mcmc

CKM_gen_rhopi:  CKM_gen_rhopi.cc
##test:  test.cc
# -----------------------------------------------------------------------------
	$(CXX) $(CXXFLAGS) -c -I. $<
	$(LD) $(LDFLAGS) -o CKM_gen_rhopi CKM_gen_rhopi.o  $(GLIBS)

CKM_gen_rhopi_mcmc:  CKM_gen_rhopi_mcmc.cc
# -----------------------------------------------------------------------------
	$(CXX) $(CXXFLAGS) -c -I. $<
	$(LD) $(LDFLAGS) -o CKM_gen_rhopi_mcmc CKM_gen_rhopi_mcmc.o  $(GLIBS)


# ================================================================================
clean:
	rm -f CKM_gen_rhopi.o  CKM_gen_rhopi
	rm -f CKM_gen_rhopi_mcmc.o  CKM_gen_rhopi_mcmc
	rm -f *~
