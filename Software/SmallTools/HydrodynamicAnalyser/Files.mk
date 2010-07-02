INCS = -I$(HDF5) -I$(BASE1INC) -I$(BASE2INC)
LIBS = $(NETCDFLIBS)
SRCS = \
       ModuleHydrodynamicAnalyser.$(F) \
       MainHydrodynamicAnalyser.$(F)
FILES = $(SRCS:.$(F)=.$(S))
METAFILES = \
            Files.smk \
            Dependencies.smk
TARGET = $(HYDROANALYSE)
