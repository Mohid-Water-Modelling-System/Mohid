# Data file

This file contains the relevant information to run the MohidReservoirOptimization program. Like most of Mohid system's files, it is uses keywords, which are:

- START - initial time 

$$
AAAA MM DD HH MM SS
$$

;
- END - final time 

$$
AAAA MM DD HH MM SS
$$

;
- DT - time step in seconds;
- W2\_BATH\_FORMAT - determines if the bathymetry file is in [free form](free-form.md) or [CQUAL-W2's format](../io/formats/cqual-w2-s-format.md);
- W2\_DATA\_FORMAT - determines if the control file is in [free form](free-form.md) or [CQUAL-W2's format](../io/formats/cqual-w2-s-format.md);
- IMX - number of segments - only if W2\_DATA\_FORMAT is set to free;
- KMX - number of layers - only if W2\_DATA\_FORMAT is set to free;
- NWB - number of bathymetry files (must be 1) - only if W2\_DATA\_FORMAT is set to free;
- NBR - number of branches - only if W2\_DATA\_FORMAT is set to free;
- BE - only if W2\_DATA\_FORMAT is set to free;
- JBDN - only if W2\_DATA\_FORMAT is set to free;
- SLOPE - only if W2\_DATA\_FORMAT is set to free;
- ELBOT - only if W2\_DATA\_FORMAT is set to free;
- DHS - only if W2\_DATA\_FORMAT is set to free;
- UHS - only if W2\_DATA\_FORMAT is set to free;
- US - only if W2\_DATA\_FORMAT is set to free;
- DS - only if W2\_DATA\_FORMAT is set to free;
- BS - only if W2\_DATA\_FORMAT is set to free;
- UQB - only if W2\_DATA\_FORMAT is set to free;
- DQB - only if W2\_DATA\_FORMAT is set to free;
- RESERVOIR\_ERROR - allowed error in %, between the final reservoir volume and the acomulated volumes file.

Sample data file:

---

START  : 2006 01 01 00 00 00
END  : 2006 01 01 00 01 00
DT  : 60
!Either "fixed" or "free"
W2\_BATH\_FORMAT  : free
!Either "fixed" or "free"
W2\_DATA\_FORMAT : fixed
IMX  : 25
KMX  : 30
NWB  : 1
NBR  : 1
BE  : 1
JBDN  : 1
SLOPE  : 0.0
ELBOT  : 0.0
DHS  : 0
UHS  : 0
US  : 2
DS  : 24
BS  : 1
UQB  : 1
DQB  : 1
!allowed error %
RESERVOIR\_ERROR : 0.00025

---
