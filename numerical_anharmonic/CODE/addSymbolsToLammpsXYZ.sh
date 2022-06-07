#!/bin/bash

# Program replaces the atom indices (1, 2, ...) with the atom symbol (Cd, Se, ...) by
# mapping the indices to symbols using the index-mass mapping in the lammpsconf.par file

# Allowed atom symbols: Cd, Zn, Se, S and Te

# Author: John P. Philbin
# Date: Oct. 30th, 2018

# Constants
mCd=112.411
mZn=65.380
mSe=78.960
mS=32.065
mP=30.974
mTe=127.60
mAs=74.922
mIn=114.820

# Get the masses from the lammpsconf.par file
nAtomTypes=( "$(grep 'atom types' lammpsconf.par | awk '{print $1+1}')" )
grep -A $nAtomTypes 'Masses' lammpsconf.par | awk '{if(NR>2){print $2}}' > tmp.dat 
m=( "$(awk '{if(NR==1){print $1}}' tmp.dat)" )
masses[0]=$m
m=( "$(awk '{if(NR==2){print $1}}' tmp.dat)" )
masses[1]=$m
m=( "$(awk '{if(NR==3){print $1}}' tmp.dat)" )
masses[2]=$m
m=( "$(awk '{if(NR==4){print $1}}' tmp.dat)" )
masses[3]=$m
m=( "$(awk '{if(NR==5){print $1}}' tmp.dat)" )
masses[4]=$m
rm tmp.dat

# Map atom masses to atom symbols
for i in {0..4}; do
	varCd=( "$(awk -v a=${masses[$i]} -v b=$mCd 'BEGIN{printf "%.3f\n", a-b}')" )
	varZn=( "$(awk -v a=${masses[$i]} -v b=$mZn 'BEGIN{printf "%.3f\n", a-b}')" )
	varSe=( "$(awk -v a=${masses[$i]} -v b=$mSe 'BEGIN{printf "%.3f\n", a-b}')" )
	varS=( "$(awk -v a=${masses[$i]} -v b=$mS 'BEGIN{printf "%.3f\n", a-b}')" )
	varP=( "$(awk -v a=${masses[$i]} -v b=$mP 'BEGIN{printf "%.3f\n", a-b}')" )
	varTe=( "$(awk -v a=${masses[$i]} -v b=$mTe 'BEGIN{printf "%.3f\n", a-b}')" )
	varAs=( "$(awk -v a=${masses[$i]} -v b=$mAs 'BEGIN{printf "%.3f\n", a-b}')" )
	varIn=( "$(awk -v a=${masses[$i]} -v b=$mIn 'BEGIN{printf "%.3f\n", a-b}')" )
	if   (( $(echo "$varCd > -0.1" | bc -l) )) && (( $(echo "$varCd < 0.1" | bc -l) )); then 
		a[$i]="Cd"
	elif (( $(echo "$varZn > -0.1" | bc -l) )) && (( $(echo "$varZn < 0.1" | bc -l) )); then 
		a[$i]="Zn"
	elif (( $(echo "$varSe > -0.1" | bc -l) )) && (( $(echo "$varSe < 0.1" | bc -l) )); then
		a[$i]="Se"
	elif (( $(echo "$varS > -0.1" | bc -l) )) && (( $(echo "$varS < 0.1" | bc -l) )); then
		a[$i]="S"
	elif (( $(echo "$varP > -0.1" | bc -l) )) && (( $(echo "$varP < 0.1" | bc -l) )); then
		a[$i]="P"
	elif (( $(echo "$varTe > -0.1" | bc -l) )) && (( $(echo "$varTe < 0.1" | bc -l) )); then
		a[$i]="Te"
	elif (( $(echo "$varAs > -0.1" | bc -l) )) && (( $(echo "$varAs < 0.1" | bc -l) )); then
		a[$i]="As"
	elif (( $(echo "$varIn > -0.1" | bc -l) )) && (( $(echo "$varIn < 0.1" | bc -l) )); then
		a[$i]="In"
	fi
done

# Add atom symbols into xyz file (passed in $1) and rename it new_$1
awk -v a=${a[0]} -v b=${a[1]} -v c=${a[2]} -v d=${a[3]} -v e=${a[4]} '{if(NR==1){print $0"\n"}else if(NR>2){if($1==1){x=a}else if($1==2){x=b}else if($1==3){x=c}else if($1==4){x=d}else if($1==5){x=e}; print x" "$2" "$3" "$4;}}' $1 > new_$1
