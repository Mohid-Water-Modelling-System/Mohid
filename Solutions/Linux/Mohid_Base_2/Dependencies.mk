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

ModuleGridData.${O} : ModuleHorizontalGrid.${O}
ModuleHorizontalMap.${O} : ModuleGridData.${O}
ModuleBoxDif.${O} : ModuleHorizontalGrid.${O}
ModuleBasinGeometry.${O} : ModuleGridData.${O}
ModuleGeometry.${O} : ModuleHorizontalMap.${O}
ModuleMap.${O} : ModuleGeometry.${O}
ModuleAdvectionDiffusion.${O} : ModuleMap.${O}
ModuleFillMatrix.${O} : ModuleBoxDif.${O} \
                     ModuleGeometry.${O}
ModuleInterpolation.${O} : ModuleMap.${O}
ModuleAtmosphere.${O} : ModuleStatistic.${O} \
                     ModuleFillMatrix.${O}

#--------End of dependencies list---------

#----------------------------------------------------------------------------------------------------------
#MOHID Water Modelling System.
#Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
#----------------------------------------------------------------------------------------------------------
