#!/usr/bin/env bash

echo "How many '''keywords''' are there in MohidWater (including Base1 and Base2 libraries)?"
for i in `find | grep MOHID | grep -v Land | grep -v River | grep F90`; do cat $i | grep -re "call GetData" ; done | wc -l
#----
echo "And in MohidWater and MohidLand together?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -re "call GetData" ; done | wc -l
#----
echo "How many '''lines of code''' is there in MohidWater and MohidLand (including Base1 and Base2)?"
for i in `find | grep MOHID | grep F90`; do cat $i; done | wc -l              
#----
echo "And without empty lines?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep "\w"; done | wc -l  
#----
echo "And without comment-lines?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep "\w" | grep -v "^ *\!"; done | wc -l
#----
echo "How many parallel '''openmp''' zones is there in the code?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep "OMP PARALLEL"; done | wc -l
#----
echo "How many subroutines is there in the code?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *subroutine"; done | wc -l
#----
echo "How many '''ifs''' are there in the code?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *if ("; done | wc -l
#----
echo "How many '''calls''' are there in the code?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *call "; done | wc -l
#----
echo "How many '''types''' are there in the code?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *type "; done | wc -l
#----
echo "How many '''allocates''' are there in the code?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *allocate "; done | wc -l
#----
echo "And '''deallocates'''?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *deallocate "; done | wc -l
#----
echo "How many '''do loops''' are there in the code?"
for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *do "; done | wc -l
#----
