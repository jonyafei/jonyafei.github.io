#! /usr/bin/sh
############################################################################
####    This script is writen by Yafei Jiang, Dec 18,2019   @Jonyafei   ####
####    Find all .log files and print corresponding energy              ####
####    current dir: ./opt/  ./freq/  ./hbs/  ./solv/                   ####
####    file name: *-opt.log  *-freq.log  *-hbs.log  *-solv.log         ####
############################################################################

cd ./opt/
files=(`ls *.log`)
names=()
Eopt=()
for i in ${!files[@]}
do
names[i]=`echo ${files[$i]} | sed 's/-opt.log//'`
Eopt[i]=`grep "SCF Done" ${files[$i]} | tail -1 | awk '{printf "%.6f\n",$5}'`
done

Ehbs=()
Ezpe=()
Gcorr=()
Esolv=()
for i in ${!names[@]}
do
    cd ../hbs/
    if [ ! -f ${names[$i]}-hbs.log ];then
    Ehbs[i]=None
    else
    Ehbs[i]=`grep "SCF Done" ${names[$i]}-hbs.log | tail -1 | awk '{printf "%.6f\n",$5}'`
    fi
    
    cd ../freq/
    if [ ! -f ${names[$i]}-freq.log ];then
    Ezpe[i]=None
    Gcorr[i]=None
    else
    Ezpe[i]=`grep "Zero-point correction" ${names[$i]}-freq.log | awk '{printf "%.6f\n",$3}'`
    Gcorr[i]=`grep "Thermal correction to Gibbs Free Energy" ${names[$i]}-freq.log | awk '{printf "%.6f\n",$7}'`
    fi
    
    cd ../solv/
    if [ ! -f ${names[$i]}-solv.log ];then
    Esolv[i]=None
    else
    Esolv[i]=`grep "SCF Done" ${names[$i]}-solv.log | tail -1 | awk '{printf "%.6f\n",$5}'`
    fi
done

cd ..
for i in ${!files[@]}
do
echo ${names[$i]} ${Eopt[$i]} ${Ezpe[$i]} ${Gcorr[$i]} ${Ehbs[$i]} ${Esolv[$i]} >> data.out
done

index="names\tE0\tEzpe\tGcorr\tEhbs\tEsolv"
#sed -r "1i$index" data.out | awk '{printf "%16s%16s%16s%16s%16s%16s\n",$1,$2,$3,$4,$5,$6}' > energy.out
sed -r "1i$index" data.out | column -t -R1,2,3,4,5,6 > energy.out
rm data.out