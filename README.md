# MOHID - Water Modelling System  

MOHID is short for Modelo Hidrodinâmico which is hydrodynamic model in Portuguese. MOHID is a three-dimensional water modelling system, developed by MARETEC (Marine and Environmental Technology Research Center) at Instituto Superior Técnico (IST) which belongs to Lisbon University.

## What is this repository?
This is the MOHID Water Modelling System OFFICIAL repository

## Overview
MOHID is a modular finite volumes water-modelling system written in ANSI-Fortran95 using an Object-oriented programming philosophy, integrating diverse mathematical models and supporting graphical user interfaces that manage all the pre- and post-processing. 
MOHID allows the adoption of an integrated modelling philosophy, not only of processes (physical and biogeochemical), but also of different scales (allowing the use of nested models) and systems (estuaries and watersheds), due to the adoption of an object oriented programming philosophy.
The development of MOHID started back in 1985. Since that time a continuous development effort of new features has been maintained. Model updates and improvements were made available in a regular basis were used in the framework of many research and engineering projects.
All programs included in MOHID Water Modelling System are built on the top of one or more base libraries and the two core executables files can be found at the top of the pyramid:
* MOHID Water – Three-dimensional mathematical model to simulate surface water bodies.
* MOHID Land – Watershed mathematical model or Hydrological transport model designed to simulate drainage basin and aquifer;

Smaller utility programs are easily built on the top of the libraries, which are usually designed for pre or post-processing results of the models. This support tools are normally managed by graphical user interfaces which allow management of input data, control of program execution, and output results analysis, along with other pre- and post-processing operations.
The integration of MOHID’s different tools can be easily achieved since these tools are based on the same framework. This coupling can thus be used to study the water cycle and its associated processes in an integrated approach.

## Help, Bugs, Feedback
If you need help with MOHID, want to keep up with progress, chat with developers or ask any other questions about MOHID, you can hang out by mail: <general@mohid.com> or consult our [MOHID wiki](http://wiki.mohid.com). You can also subscribe to our [MOHID forum](http://forum.mohid.com). To report bugs, please create a GitHub issue or contact any developers. More information consult <http://www.mohid.com>

## License
GNU General Public License. See the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) web page for more information.
