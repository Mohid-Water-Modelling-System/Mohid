#------------------------------------------------------------------------------
#        IST/MARETEC, Water Modelling Group, Mohid modelling system
#------------------------------------------------------------------------------
#
# TITLE         : Mohid Model
# PROJECT       : Mohid Base 1
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

INCS = -I$(HDF5MODINC)
LIBS =

ifeq ($(IS_PROJ4F),true)
    INCS := $(INCS) -I$(PROJ4FINC)
    LIBS := $(LIBS) $(PROJ4FLIBS)
endif

SRCDIR = $(SRCREP)/MOHIDBase1
FILES = \
        ModuleGlobalData.$(S) \
        ModuleTime.$(S) \
        ModuleEnterData.$(S) \
        ModuleFunctions.$(S) \
	ModuleBivalve.$(S) \
	ModuleSeagrassSedimInteraction.$(S) \
	ModuleSeagrassWaterInteraction.$(S) \
        ModuleBenthos.$(S) \
        ModuleBenthicEcology.$(S) \
        ModuleCEQUALW2.$(S) \
        ModuleMacroAlgae.$(S) \
        ModuleTimeSerie.$(S) \
        ModuleDischarges.$(S) \
        ModuleDrawing.$(S) \
        ModuleLUD.$(S) \
        ModuleWaterQuality.$(S) \
        ModuleSedimentQuality.$(S) \
        ModuleLife.$(S) \
        ModuleInterface.$(S) \
        ModuleHydroIntegration.$(S) \
        ModuleLightExtinction.$(S) \
        ModuleStopWatch.$(S) \
        ModuleTriangulation.$(S) \
        ModuleHDF5.$(S) \
        ModuleDrainageNetwork.$(S) \
        ModuleProfile.$(S) \
        ModuleTask2000.F \
        ModuleWWTPQ.$(S)
#	ModuleCuda.$(S) 

ifeq ($(IS_PROJ4F),true)
    FILES := proj4.${S} ${FILES} 
endif

TARGET = $(BASE1)

#--------Makefile rules----------------

SHELL = /bin/sh

include ../Makefiles/Makefiles.mk

#--------Dependencies list-------------

ModuleTime.${Obj} 			: ModuleGlobalData.${Obj}
ModuleStopWatch.${Obj} 			: ModuleTime.${Obj}
ModuleEnterData.${Obj} 			: ModuleStopWatch.${Obj}
ModuleFunctions.${Obj} 			: ModuleEnterData.${Obj}
ModuleSeagrassSedimInteraction.${Obj} 	: ModuleFunctions.${Obj} 
ModuleSeagrassWaterInteraction.${Obj} 	: ModuleFunctions.${Obj} 
ModuleBivalve.${Obj} 			: ModuleEnterData.${Obj}
ModuleBenthos.${Obj} 			: ModuleEnterData.${Obj}
ModuleCEQUALW2.${Obj} 			: ModuleFunctions.${Obj}
ModuleMacroAlgae.${Obj} 		: ModuleFunctions.${Obj}
ModuleTimeSerie.${Obj} 			: ModuleEnterData.${Obj}
ModuleDischarges.${Obj} 		: ModuleFunctions.${Obj} \
                    			  ModuleTimeSerie.${Obj} \
                    			  ModuleDrawing.${Obj}
ModuleDrawing.${Obj} 			: ModuleFunctions.${Obj}
ModuleLUD.${Obj} 			: ModuleGlobalData.${Obj}
ModuleWaterQuality.${Obj} 		: ModuleFunctions.${Obj} \
                      			  ModuleLUD.${Obj}
ModuleSedimentQuality.${Obj} 		: ModuleFunctions.${Obj} \
                        		  ModuleLUD.${Obj}
ModuleBenthicEcology.${Obj} 		: ModuleFunctions.${Obj}
ModuleWWTPQ.${Obj} 			: ModuleFunctions.${Obj} \
					  ModuleLUD.${Obj}	
ModuleLife.${Obj} 			: ModuleFunctions.${Obj}
ModuleInterface.${Obj} 			: ModuleWaterQuality.${Obj} \
					  ModuleSedimentQuality.${Obj} \
					  ModuleCEQUALW2.${Obj} \
					  ModuleMacroAlgae.${Obj} \
					  ModuleLife.${Obj} \
					  ModuleBenthos.${Obj} \
					  ModuleBenthicEcology.${Obj} \
					  ModuleWWTPQ.${Obj} \
					  ModuleSeagrassSedimInteraction.${Obj} \
					  ModuleSeagrassWaterInteraction.${Obj} \
					  ModuleBivalve.${Obj}
ModuleHydroIntegration.${Obj} 		: ModuleTime.${Obj}
ModuleLightExtinction.${Obj}		: ModuleFunctions.${Obj} \
                          		  ModuleTimeSerie.${Obj}
ModuleTriangulation.${Obj} 		: ModuleGlobalData.${Obj}
ModuleHDF5.${Obj} 			: ModuleGlobalData.${Obj}
ModuleDrainageNetwork.${Obj}	 	: ModuleDischarges.${Obj} \
                          		  ModuleLightExtinction.${Obj} \
                          		  ModuleStopWatch.${Obj} \
                          		  ModuleInterface.${Obj}
#ModuleCuda.${Obj} 			: ModuleGlobalData.${Obj}
ModuleProfile.${Obj} 			: ModuleEnterData.${Obj} \
                  			  ModuleHDF5.${Obj}

ModuleTask2000.${Obj} :
	@$(CC) $(CCFLAGS) $(INCS) $(SRCDIR)/ModuleTask2000.F
	@echo $* .................. [OK]    


#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------


