#------------------------------------------------------------
#
# Please define an environnement variable PATHTUTOMIP
# equal to the XXX directory where Lemon,Graphviz and Cplex have been installed
# In order to these:
# edit ~/.bashrc add "export PATHTUTOMIP=XXX"
# close your xterm and star new ones
#
#------------------------------------------------------------

CPLEXDIR      = $(PATHTUTOMIP)/cplex-12.7.1/cplex
CONCERTDIR    = $(PATHTUTOMIP)/cplex-12.7.1/concert
LEMONDIR      = $(PATHTUTOMIP)/lemon-1.3.1


SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic


# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0


# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wno-ignored-attributes

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LEMONLIBDIR	= $(LEMONDIR)/lib

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -L$(LEMONLIBDIR)


CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
LEMONINCDIR   = $(LEMONDIR)/include

CCLNFLAGS = -lconcert -lilocplex -lcplex -lemon -lm -lpthread

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(LEMONINCDIR)




# ---------------------------------------------------------------------
# Your part
# ---------------------------------------------------------------------



all: BC_Acyclic

Graph.o: ../../Graph/Graph.cpp ../../Graph/Graph.h
	$(CCC) -c $(CCFLAGS) $(CCLNDIRS) -c ../../Graph/Graph.cpp -o ../../Graph/Graph.o

CircuitSeparation.o: CircuitSeparation.cpp ../../Graph/Graph.h
	$(CCC) -c $(CCFLAGS) $(CCLNDIRS) CircuitSeparation.cpp

BC_Acyclic.o: BC_Acyclic.cpp
	$(CCC) -c $(CCFLAGS) $(CCLNDIRS) BC_Acyclic.cpp


BC_Acyclic: BC_Acyclic.o CircuitSeparation.o ../../Graph/Graph.o 
	$(CCC) $(CCFLAGS) $(CCLNDIRS) BC_Acyclic.o CircuitSeparation.o ../../Graph/Graph.o   -o ../../bin/BC_Acyclic $(CCLNFLAGS)


clean:
	rm -f *.o ../../bin/BC_Acyclic
