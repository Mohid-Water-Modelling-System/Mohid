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

ModuleGauge.${Obj} : ModuleToga.${Obj}
ModuleSand.${Obj} : ModuleWaves.${Obj}
ModuleSedimentProperties.${Obj} : ModuleConsolidation.${Obj}
ModuleOpenBoundary.${Obj} : ModuleGauge.${Obj}
ModuleTurbGOTM.${Obj} : ModuleGOTM.${Obj}
ModuleTurbulence.${Obj} : ModuleTurbGOTM.${Obj}
ModuleHydrodynamic.${Obj} : ModuleAssimilation.${Obj} \
                       ModuleHydrodynamicFile.${Obj} \
                       ModuleOpenBoundary.${Obj} \
                       ModuleTurbulence.${Obj}
ModuleWaterProperties.${Obj} : ModuleFreeVerticalMovement.${Obj} \
                          ModuleHydrodynamic.${Obj}
ModuleLagrangian.${Obj} : ModuleJet.${Obj} \
                     ModuleOil.${Obj} \
                     ModuleOil_0D.${Obj} \
                     ModuleWaves.${Obj} \
                     ModuleWaterProperties.${Obj}
ModuleLagrangianGlobal.${Obj} : ModuleJet.${Obj} \
                     ModuleOil.${Obj} \
                     ModuleOil_0D.${Obj} \
                     ModuleWaves.${Obj} \
                     ModuleWaterProperties.${Obj}                    	
ModuleInterfaceWaterAir.${Obj} : ModuleLagrangian.${Obj} \
                               ModuleLagrangianGlobal.${Obj} \

ModuleInterfaceSedimentWater.${Obj} : ModuleSand.${Obj} \
                                 ModuleLagrangian.${Obj} \
                                 ModuleLagrangianGlobal.${Obj} \
                                 ModuleSedimentProperties.${Obj}
ModuleModel.${Obj} : ModuleInterfaceWaterAir.${Obj} \
                ModuleInterfaceSedimentWater.${Obj}
ModuleMain.${Obj} : ModuleModel.${Obj}

#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------
