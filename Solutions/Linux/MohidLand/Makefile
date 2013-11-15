#------------------------------------------------------------------------------
#        IST/MARETEC, Water Modelling Group, Mohid modelling system
#------------------------------------------------------------------------------
#
# TITLE         : Mohid Model
# PROJECT       : Mohid Land
# URL           : http://www.mohid.com
# AFFILIATION   : IST/MARETEC, Marine Modelling Group
# DATE          : November 2013
# REVISION      : Eduardo Jauch - v1.1
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

INCS = -I$(HDF5MODINC) -I$(BASE1INC) -I$(BASE2INC)
LIBS = $(BASELIBS) $(HDFLIBS)
SRCDIR = $(SRCREP)/MOHIDLand
FILES = \
    ModuleBasin.$(S) \
    ModulePorousMedia.$(S) \
    ModulePorousMediaProperties.$(S) \
    ModuleRunOff.$(S) \
    ModuleRunoffProperties.$(S) \
    ModuleVegetation.$(S) \
    MohidLand.$(S)
TARGET = MohidLand$(SUFFPROG)

#----------Makefile rules--------------

SHELL = /bin/sh 

include ../Makefiles/Makefiles.mk

#--------Dependencies list-------------

ModulePorousMediaProperties.${Obj} 	: ModulePorousMedia.${Obj}
ModuleRunoffProperties.${Obj} 		: ModuleRunOff.${Obj}
ModuleBasin.${Obj} 			: ModuleRunoffProperties.${Obj} \
                    			  ModulePorousMediaProperties.${Obj} \
                    			  ModuleVegetation.${Obj}

#--------End of dependencies list---------

#---------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#---------------------------------------------------------------------------------------------------
