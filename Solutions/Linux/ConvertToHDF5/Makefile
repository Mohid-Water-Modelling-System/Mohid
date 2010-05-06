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

INCS = -I$(HDF5) -I$(BASE1INC) -I$(BASE2INC)
LIBS = $(BASELIBS)
SRCDIR = $(SRCCONVERT)
FILES := 
FILES := $(FILES) ConvertToHDF5.$(S)
#FILES := $(FILES) ModuleAladinFormat.$(S)
#FILES := $(FILES) ModuleARPSFormat.$(S)
#FILES := $(FILES) ModuleARPSToWW3.$(S)
#FILES := $(FILES) ModuleCFFormat.$(S)
#FILES := $(FILES) ModuleCFPolcomormat.$(S)
#FILES := $(FILES) ModuleConvertModisL2.$(S)
#FILES := $(FILES) ModuleConvertModisL3Mapped.$(S)
#FILES := $(FILES) ModuleConvertOceanColorL2.$(S)
#FILES := $(FILES) ModuleCOWAMAASCIIWind.$(S)
#FILES := $(FILES) ModuleERA40Format.$(S)
FILES := $(FILES) ModuleEUCenterFormat.$(S)
FILES := $(FILES) ModuleGFSasciiWind.$(S)
FILES := $(FILES) ModuleGlueHDF5Files.$(S)
FILES := $(FILES) ModuleHDF5ToASCIIandBIN.$(S)
FILES := $(FILES) ModuleHellermanRosensteinAscii.$(S)
#FILES := $(FILES) ModuleHYCOMFormat.$(S)
FILES := $(FILES) ModuleInterpolateGrids.$(S)
FILES := $(FILES) ModuleLevitusFormat.$(S)
#FILES := $(FILES) ModuleMERCATORFormat.$(S)
FILES := $(FILES) ModuleMM5Format.$(S)
#FILES := $(FILES) ModuleMOG2DFormat.$(S)
FILES := $(FILES) ModulePatchHDF5Files.$(S)
FILES := $(FILES) ModuleSWAN.$(S)
FILES := $(FILES) ModuleTecnoceanAscii.$(S)
#FILES := $(FILES) ModuleWOAFormat.$(S)
#FILES := $(FILES) ModuleWRFFormat.$(S)
#FILES := $(FILES) ncdflib.$(S)

TARGET = $(CONVERT)

#----------Makefile rules--------------

SHELL = /bin/sh 

include ../Makefiles/Makefiles.mk

#--------Dependencies list-------------

CONVERTDEPEND := 
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleAladinFormat.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleARPSFormat.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleARPSToWW3.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleCFFormat.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleCFPolcomormat.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleConvertModisL2.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleConvertModisL3Mapped.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleConvertOceanColorL2.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleCOWAMAASCIIWind.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleERA40Format.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleEUCenterFormat.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleGFSasciiWind.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleGlueHDF5Files.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleHDF5ToASCIIandBIN.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleHellermanRosensteinAscii.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleHYCOMFormat.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleInterpolateGrids.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleLevitusFormat.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleMERCATORFormat.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleMM5Format.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleMOG2DFormat.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModulePatchHDF5Files.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleSWAN.$(O)
CONVERTDEPEND := $(CONVERTDEPEND) ModuleTecnoceanAscii.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleWOAFormat.$(O)
#CONVERTDEPEND := $(CONVERTDEPEND) ModuleWRFFormat.$(O)

#ModuleCFFormat.${O} : ncdflib.${O}
ConvertToHDF5.${O} : $(CONVERTDEPEND)

#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------
