#------------------------------------------------------------------------------
#        IST/MARETEC, Water Modelling Group, Mohid modelling system
#------------------------------------------------------------------------------
#
# TITLE         : Mohid Model
# PROJECT       : Mohid Water
# URL           : http://www.mohid.com
# AFFILIATION   : IST/MARETEC, Marine Modelling Group
# DATE          : September 2006
# REVISION      : Guillaume Riflet - v1.0
# DESCRIPTION   : Makefile
#------------------------------------------------------------------------------
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License 
#version 2, as published by the Free Software Foundation.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#------------------------------------------------------------------------------

INCS := -I$(HDF5INC) -I$(BASE1INC) -I$(BASE2INC)
LIBS = $(BASELIBS) $(HDFLIBS)

ifeq ($(IS_NETCDF),true)
    INCS := $(INCS) -I$(NETCDFINC)
    LIBS := $(LIBS) $(NETCDFLIBS)
endif

ifeq ($(IS_PROJ4F),true)
    INCS := $(INCS) -I$(PROJ4FINC)
    LIBS := $(LIBS) $(PROJ4FLIBS)
endif

SRCDIR = $(SRCREP)/ConvertToHDF5

NCFILES :=
NCFILES := $(NCFILES) ModuleAladinFormat.$(S)
NCFILES := $(NCFILES) ModuleARPSFormat.$(S) 
NCFILES := $(NCFILES) ModuleARPSToWW3.$(S)
NCFILES := $(NCFILES) ModuleCFFormat.$(S)
NCFILES := $(NCFILES) ModuleCFPolcomFormat.$(S)
NCFILES := $(NCFILES) ModuleCOWAMAAsciiWind.$(S)
NCFILES := $(NCFILES) ModuleERA40Format.$(S)
NCFILES := $(NCFILES) ModuleHYCOMFormat.$(S)
NCFILES := $(NCFILES) ModuleMERCATORFormat.$(S)
NCFILES := $(NCFILES) ModuleMOG2DFormat.$(S) 
NCFILES := $(NCFILES) ModuleWOAFormat.$(S)
NCFILES := $(NCFILES) ncdflib.$(S)

PROJ4FILES :=
PROJ4FILES := $(PROJ4FILES) ModuleWRFFormat.$(S)

FILES := 
FILES := $(FILES) ConvertToHDF5.$(S)
#FILES := $(FILES) ModuleConvertModisL2.$(S)
#FILES := $(FILES) ModuleConvertModisL3Mapped.$(S)
#FILES := $(FILES) ModuleConvertOceanColorL2.$(S)
FILES := $(FILES) ModuleEUCenterFormat.$(S)
FILES := $(FILES) ModuleGFSasciiWind.$(S)
FILES := $(FILES) ModuleGlueHDF5Files.$(S)
FILES := $(FILES) ModuleHDF5ToASCIIandBIN.$(S)
FILES := $(FILES) ModuleHellermanRosensteinAscii.$(S)
FILES := $(FILES) ModuleInterpolateGrids.$(S)
FILES := $(FILES) ModuleLevitusFormat.$(S)
FILES := $(FILES) ModuleMM5Format.$(S)
FILES := $(FILES) ModulePatchHDF5Files.$(S)
FILES := $(FILES) ModuleSWAN.$(S)
FILES := $(FILES) ModuleTecnoceanAscii.$(S)

ifeq ($(IS_NETCDF),true)
    FILES := $(FILES) $(NCFILES)
endif

ifeq ($(IS_PROJ4F),true)
   FILES := $(FILES) $(PROJ4FILES)
endif

TARGET = ConvertToHDF5$(SUFFPROG)

#----------Makefile rules--------------

SHELL = /bin/sh 

include ../Makefiles/Makefiles.mk

#--------Dependencies list-------------

NCCONVERTDEPEND :=
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleAladinFormat.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleARPSFormat.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleARPSToWW3.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleCFFormat.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleCFPolcomFormat.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleCOWAMAAsciiWind.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleERA40Format.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleHYCOMFormat.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleMERCATORFormat.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleMOG2DFormat.$(O)
NCCONVERTDEPEND := $(NCCONVERTDEPEND) ModuleWOAFormat.$(O)

CONVERTDEPEND := 
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleConvertModisL2.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleConvertModisL3Mapped.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleConvertOceanColorL2.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleEUCenterFormat.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleGFSasciiWind.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleGlueHDF5Files.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleHDF5ToASCIIandBIN.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleHellermanRosensteinAscii.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleInterpolateGrids.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleLevitusFormat.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleMM5Format.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModulePatchHDF5Files.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleSWAN.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleTecnoceanAscii.$(O)

PROJ4CONVERTDEPEND :=
PROJ4CONVERTDEPEND := $(PROJ4CONVERTDEPEND) ModuleWRFFormat.$(O)

ifeq ($(IS_NETCDF),true)
    CONVERTDEPEND := $(CONVERTDEPEND) $(NCCONVERTDEPEND)
endif

ifeq ($(IS_PROJ4F),true)
    CONVERTDEPEND := $(CONVERTDEPEND) $(PROJ4CONVERTDEPEND)
endif

## ---dependencies actually are only defined here -----
ifeq ($(IS_NETCDF),true)
    $(NCCONVERTDEPEND) : ncdflib.$(O)
endif
ConvertToHDF5.${O} : $(CONVERTDEPEND)
## --- end of actual dependencies

#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------
