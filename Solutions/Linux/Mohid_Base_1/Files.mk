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

INCS = -I$(HDF5)
LIBS =
SRCDIR = $(SRCBASE1)
FILES = \
        ModuleGlobalData.$(S) \
        ModuleTime.$(S) \
        ModuleEnterData.$(S) \
        ModuleFunctions.$(S) \
        ModuleBenthos.$(S) \
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
        ModuleProfile.$(S)
METAFILES = \
        Files.smk \
        Dependencies.smk
TARGET = $(BASE1)

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------
