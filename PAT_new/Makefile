APP_NAME := anapp120pions
SOURCE_FILES :=  heditor.cc hcommondef.cc hcut.cc htrackcut.cc htimecut.cc hgraphcut.cc hdedxcut.cc hparticlecandidate.cc hparticle.cc hhypcandidate.cc heventpool.cc hpool.cc hparticledatapool.cc hparticlepool.cc hhypdatapool.cc hhyppool.cc hpidpool.cc htrackplayer.cc hparticleplayer.cc hhypplayer.cc hntuple.cc hfile.cc houtputfile.cc houtput.cc main.cc wall/fwpad.cc wall/fwdetector.cc wall/fwsinglehit.cc wall/fwclusterhit.cc wall/fwhit.cc

USES_RFIO    := no
USES_CERNLIB := no
USES_ORACLE  := no
#CPP_FLAGS    += -g -fno-inline -DDEBUG
CPP_FLAGS    +=  # -g # -m32 


include $(HADDIR)/hades.def.mk

#LD	     += -g -m32 -fno-inline -DDEBUG
LD	     += # -g # -m32 

# override default list of linked Hydra libraries - before they can act on the rules
#HYDRA_LIBS += -lHyp
HYDRA_LIBS := \
        -lHyp \
        -lRich \
        -lRichUtil \
        -lMdcTrackD \
        -lMdc \
        -lMdcUtil \
        -lMdcTrackS \
        -lMdcTrackG \
        -lKick \
        -lShower \
        -lShowerUtil \
        -lShowerTofino \
        -lStart \
        -lTof \
        -lTofino \
        -lTrigger \
        -lTriggerUtil \
        -lPid \
        -lPidUtil \
        -lPairs \
        -lWall \
        -lHydra

.PHONY:  default
default: build install

include $(HADDIR)/hades.app.mk
