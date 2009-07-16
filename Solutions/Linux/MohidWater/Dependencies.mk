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
# DESCRIPTION   : Makefile dependencies list
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
                     ModuleWaves.${O} \
                     ModuleWaterProperties.${O}
ModuleInterfaceWaterAir.${O} : ModuleLagrangian.${O}
ModuleInterfaceSedimentWater.${O} : ModuleSand.${O} \
                                 ModuleLagrangian.${O} \
                                 ModuleSedimentProperties.${O}
ModuleModel.${O} : ModuleInterfaceWaterAir.${O} \
                ModuleInterfaceSedimentWater.${O}
ModuleMain.${O} : ModuleModel.${O}

#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------
