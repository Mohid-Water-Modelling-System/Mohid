INCS = -I$(HDF5) -I$(BASE1INC) -I$(BASE2INC) -I$(NETCDFINC)
LIBS = $(NETCDFLIBS)
SRCS = \
       ModuleNETCDF.$(F) \
       ConvertToNETCDF.$(F)
FILES = $(SRCS:.$(F)=.$(S))
METAFILES = \
            Files.smk \
            Dependencies.smk
TARGET = $(CONVERTTONETCDF)
SOSPROJ = SmallTools/$(basename $(TARGET))
