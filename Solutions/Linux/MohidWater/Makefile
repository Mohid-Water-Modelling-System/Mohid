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

INCS = -I$(HDF5INC) -I$(BASE1INC) -I$(BASE2INC)
LIBS = $(BASELIBS) $(HDFLIBS)
SRCDIR = $(SRCREP)/MOHIDWater
FILES = \
       ModuleAssimilation.$(S) \
       ModuleConsolidation.$(S) \
       ModuleFreeVerticalMovement.$(S) \
       ModuleJet.$(S) \
       ModuleHydrodynamicFile.$(S) \
       ModuleOil.$(S) \
       ModuleOil_0D.$(S) \
       ModuleToga.$(S) \
       ModuleWaves.$(S) \
       ModuleGauge.$(S) \
       ModuleSand.$(S) \
       ModuleSedimentProperties.$(S) \
       ModuleOpenBoundary.$(S) \
       ModuleGOTM.$(S) \
       ModuleTurbGOTM.$(S) \
       ModuleTurbulence.$(S) \
       ModuleHydrodynamic.$(S) \
       ModuleWaterProperties.$(S) \
       ModuleLagrangian.$(S) \
       ModuleLagrangianGlobal.$(S) \
       ModuleInterfaceWaterAir.$(S) \
       ModuleInterfaceSedimentWater.$(S) \
       ModuleModel.$(S) \
       Main.$(S)
TARGET = MohidWater$(SUFFPROG)

#----------Makefile rules--------------

SHELL = /bin/sh 

include ../Makefiles/Makefiles.mk

#--------Dependencies list-------------

ModuleGauge.${O} : ModuleToga.${O}
ModuleSand.${O} : ModuleWaves.${O}
ModuleSedimentProperties.${O} : ModuleConsolidation.${O}
ModuleOpenBoundary.${O} : ModuleGauge.${O}
ModuleTurbGOTM.${O} : ModuleGOTM.${O}
ModuleTurbulence.${O} : ModuleTurbGOTM.${O}
ModuleHydrodynamic.${O} : ModuleAssimilation.${O} \
                       ModuleHydrodynamicFile.${O} \
                       ModuleOpenBoundary.${O} \
                       ModuleTurbulence.${O}
ModuleWaterProperties.${O} : ModuleFreeVerticalMovement.${O} \
                          ModuleHydrodynamic.${O}
ModuleLagrangian.${O} : ModuleJet.${O} \
                     ModuleOil.${O} \
                     ModuleOil_0D.${O} \
                     ModuleWaves.${O} \
                     ModuleWaterProperties.${O}
ModuleLagrangianGlobal.${O} : ModuleJet.${O} \
                     ModuleOil.${O} \
                     ModuleOil_0D.${O} \
                     ModuleWaves.${O} \
                     ModuleWaterProperties.${O}                    	
ModuleInterfaceWaterAir.${O} : ModuleLagrangian.${O} \
                               ModuleLagrangianGlobal.${O} \

ModuleInterfaceSedimentWater.${O} : ModuleSand.${O} \
                                 ModuleLagrangian.${O} \
                                 ModuleLagrangianGlobal.${O} \
                                 ModuleSedimentProperties.${O}
ModuleModel.${O} : ModuleInterfaceWaterAir.${O} \
                ModuleInterfaceSedimentWater.${O}
ModuleMain.${O} : ModuleModel.${O}

#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------
