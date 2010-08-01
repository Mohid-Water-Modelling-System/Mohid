#------------------------------------------------------------------------------
#        IST/MARETEC, Water Modelling Group, Mohid modelling system
#------------------------------------------------------------------------------
#
# TITLE         : Mohid Model
# PROJECT       : Mohid Water
# URL           : http://www.mohid.com
# AFFILIATION   : IST/MARETEC, Marine Modelling Group
# DATE          : September 2006
# REVISION      : Ana Rosa Trancoso
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
LIBS = $(BASELIBS) $(HDFLIBS)

ifeq ($(IS_NETCDF),true)
     INCS := $(INCS) -I$(NETCDFINC)
     LIBS := $(LIBS) $(NETCDFLIBS)
endif

SRCDIR = $(SRCREP)/ConvertToXYZ
FILES := \
    ModuleASCII.$(S) \
    ModuleEtopo2.$(S) \
    ModuleEtopo5.$(S) \
    ModuleNASA.$(S) \
    ModuleNOAAShoreLine.$(S) \
    ModuleSRTM30.$(S) \
    ConvertToXYZ.$(S)
NCFILES := \
    ModuleGEBCO.$(S)

ifeq ($(IS_NETCDF),true)
   FILES := $(FILES) $(NCFILES)
endif

TARGET = ConvertToXYZ$(SUFFPROG)

#----------Makefile rules--------------

SHELL = /bin/sh 

include ../Makefiles/Makefiles.mk

#--------Dependencies list-------------

CONVERTDEPEND := ModuleASCII.$(O) \
       		ModuleEtopo2.$(O) \
		ModuleEtopo5.$(O) \
		ModuleNASA.$(O) \
		ModuleNOAAShoreLine.$(O) \
		ModuleSRTM30.$(O)
NCCONVERTDEPEND := \
		ModuleGEBCO.$(O)

ifeq ($(IS_NETCDF),true)
    CONVERTDEPEND := $(CONVERTDEPEND) $(NCCONVERTDEPEND)
endif

ConvertToXYZ.$(O) : $(CONVERTDEPEND)

#--------End of dependencies list---------

#---------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#---------------------------------------------------------------------------------------------------
