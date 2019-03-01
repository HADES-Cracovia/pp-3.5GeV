APP_NAME := $(NAZWA)

ifdef PLIK
CPP_FLAGS+= -DPLIK=$(PLIK)
endif

SOURCE_FILES := analysisElastic.cc
USES_RFIO := no
USES_ORACLE := yes
CPP_FLAGS+= -g -I/lustre/nyx/hades/user/przygoda/PION/HNTUPLE64 

include $(HADDIR)/hades.def.mk

LD+= -g `root-config  --libs` -L/lustre/nyx/hades/user/przygoda/PION/HNTUPLE64 -lTuple 

.PHONY:  default
default: build install

include $(HADDIR)/hades.app.mk


