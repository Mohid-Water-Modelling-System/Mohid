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

#--------Dependencies list-------------

ModuleTime.${O} : ModuleGlobalData.${O}
ModuleEnterData.${O} : ModuleTime.${O}
ModuleFunctions.${O} : ModuleEnterData.${O}
ModuleBenthos.${O} : ModuleEnterData.${O}
ModuleCEQUALW2.${O} : ModuleFunctions.${O}
ModuleMacroAlgae.${O} : ModuleFunctions.${O}
ModuleTimeSerie.${O} : ModuleEnterData.${O}
ModuleDischarges.${O} : ModuleFunctions.${O} \
                     ModuleTimeSerie.${O} \
                     ModuleDrawing.${O}
ModuleDrawing.${O} : ModuleFunctions.${O}
ModuleLUD.${O} : ModuleGlobalData.${O}
ModuleWaterQuality.${O} : ModuleFunctions.${O} \
                       ModuleLUD.${O}
ModuleSedimentQuality.${O} : ModuleFunctions.${O} \
                          ModuleLUD.${O}
ModuleLife.${O} : ModuleFunctions.${O}
ModuleInterface.${O} : ModuleWaterQuality.${O} \
                    ModuleSedimentQuality.${O} \
                    ModuleCEQUALW2.${O} \
                    ModuleMacroAlgae.${O} \
                    ModuleLife.${O} \
                    ModuleBenthos.${O}
ModuleHydroIntegration.${O} : ModuleTime.${O}
ModuleLightExtinction.${O} : ModuleFunctions.${O} \
                          ModuleTimeSerie.${O}
ModuleStopWatch.${O} : ModuleTime.${O}
ModuleTriangulation.${O} : ModuleGlobalData.${O}
ModuleHDF5.${O} : ModuleGlobalData.${O}
ModuleDrainageNetwork.${O} : ModuleDischarges.${O} \
                          ModuleLightExtinction.${O} \
                          ModuleStopWatch.${O} \
                          ModuleInterface.${O}
ModuleProfile.${O} : ModuleEnterData.${O} \
                  ModuleHDF5.${O}

#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------


