# How to create a Virtual Machine from a copy

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

1. Copy .vmx (1 config) file and .vmdk (2 disk files) from the origin path to the VirtualMachines datastore.
2. In VMWare Server, on Virtual Machines Menu click "Add Virtual Machine To Inventory" choosing the .vmx file copied.(a)
3. Power on the virtual machine
4. Will ask if machine was copied. Accept this option.
5. Change machine name (and ip) inside the console in VMWARE. Restart the machine in the console so that changes take effect.
6. Machine is available for remote connections.

(a) How to Use VMWare:

1. To install VMWare server go to \\\\kepler\\Software\\VMWare\\Server.
2. If the program does not run (happened in one machine) try to go to AdministrativeTools->Services and start VMWare Host Agent
3. Initiate the program (it goes to internet browser) and sellect proceed (not recommended option in the browser) and accept the certificates. If the browser complains about scripts add the site to the trusted sites in Tools->InternetOptions->Security.
4. Give administrator credentials.
