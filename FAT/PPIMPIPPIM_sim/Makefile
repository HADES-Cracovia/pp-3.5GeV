CC=g++
TARGET=analpp
OBJS=main.o data.o PPimPipPim.o PPim.o
HEADERS=data.h PPimPipPim.h PPim.h
ROOTLIBS=`root-config  --libs` -lTuple -lHFilter -lParticle -lHydra
ROOTLIBS += -lPid -lPidUtil -lPairs -lParticle
ROOTLIBS += -lMdcTrackG -lMdcTrackD -lMdcTrackS -lKick
#HYDRA_LIBS += -lHodo -lWall -lShowerUtil -lShowerTofino -lShower -lTofino -lTof
ROOTLIBS += -lHodo -lWall -lShowerTofino -lShower -lTofino -lTof
ROOTLIBS += -lMdc -lRichUtil -lRich -lRpc -lHodo
ROOTLIBS += -lTriggerUtil -lTrigger -lStart
ROOTLIBS += -lQA -lTools
ROOTLIBS += -lTMVA





CPPFLAGS+=-g -L/u/knowakow/HNTUPLE64 -L/u/knowakow/HFILTER64 -L/u/knowakow/PLUTO64/ -L/lustre/nyx/hades/user/kempter/svn/hydra_BT/lib  -Wl,--no-as-needed
FFLAGS+=-g -I/u/knowakow/HNTUPLE64 -I/u/knowakow/HFILTER64 -I/u/knowakow/PLUTO64/src -I/lustre/nyx/hades/user/kempter/svn/hydra_BT/include
#FFLAGS+=-g -pg -fprofile-arcs
#CPPFLAGS+=-g 

ifdef NOCUT
  FFLAGS += -DNOCUT
endif

ifdef RECTANG
  FFLAGS += -DRECTANG
endif

ifdef TARG
  FFLAGS += -DTARG
endif



.SUFFIXES   : .o .cc
.SUFFIXES   : .o .C

.cc.o :
	$(CC) $(FFLAGS) `root-config  --cflags` -c $<
.C.o :
	$(CC) $(FFLAGS) `root-config  --cflags` -c $<


all: $(OBJS) $(HEADERS)
	$(CC) $(OBJS) $(CPPFLAGS) $(ROOTLIBS) -o $(TARGET) 

clean:
	-rm -rf *.o *.d $(TARGET) 

