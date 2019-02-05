copy ..\data\nomfich_1.dat ..\exe\nomfich.dat
mpiexec -np  4 MohidWater_mpi > display1.txt
DomainConsolidation.exe
copy ..\data\nomfich_2.dat ..\exe\nomfich.dat
mpiexec -np  4 MohidWater_mpi > display2.txt
DomainConsolidation.exe
pause