ARCH         := $(shell root-config --arch)
PLATFORM      = $(ARCH)

ObjSuf        = o
SrcSuf        = cc
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
USES_RFIO    := no
USES_CERNLIB := no
USES_ORACLE  := no


HYDRALIBS     = 

HYDRAINC = -I$(HADDIR)/include

PLUTOLIBS	= -L$(PLUTO_DIR) -lPluto

PLUTOINC =  -I$(PLUTO_DIR)/src -I$(PLUTO_DIR)/plugins/dalitz_mod

CXX           = g++
CXXFLAGS      = -O -Wall -g -fPIC  $(HYDRAINC) $(PLUTOINC)
LD            = g++
LDFLAGS       = -O 
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------


O1      = pluto_gen.$(ObjSuf)
S1      = pluto_gen.$(SrcSuf)

E       = pluto_gen_E

OBJS          = $(O1)
PROGRAMS      = $(E)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

clean:      
	@rm -f $(OBJS) core

$(E):       $(O1)
	$(LD) $(LDFLAGS) $^ $(HYDRALIBS) $(PLUTOLIBS) $(LIBS) $(GLIBS) $(OutPutOpt)$@
	@echo "$@ done"


.SUFFIXES: .$(SrcSuf)

#------------------------------------------------------------------------------

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
