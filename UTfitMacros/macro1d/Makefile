ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
make1Dplots: make1Dplots.cc
# ------------------------------------------------------------------------------------
	g++ -g ${ROOTCFLAGS} $@.cc  ${ROOTGLIBS} -o $@
# ====================================================================================
makepull: makepull.cc
# ------------------------------------------------------------------------------------
	g++ -g ${ROOTCFLAGS} $@.cc  ${ROOTGLIBS} -o $@
# ====================================================================================
clean:
	rm -f *~
	rm -f make1Dplots.o make1Dplots
