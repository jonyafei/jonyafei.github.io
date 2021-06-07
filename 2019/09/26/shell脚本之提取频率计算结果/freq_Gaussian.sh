#!/bin/sh

for dir in `find ./ -type d | grep 'freq'`; do
if [ ! -f $dir/*.log ];then
continue
else

echo "file in $dir is processing" 
echo $dir | sed 's/freq//' | sed 's/^..//' | awk '{printf "%-30s\n",$1}' >> name
cat $dir/*.log | grep "SCF Done" | awk '{printf "%.6f\n",$5}' >> E0
cat $dir/*.log | grep "Zero-point correction" | awk '{printf "%.6f\n",$3}' >> Ezpe
cat $dir/*.log | grep "Thermal correction to Enthalpy" | awk '{printf "%.6f\n",$5}' >> Hcorr
cat $dir/*.log | grep "Thermal correction to Gibbs Free Energy" | awk '{printf "%.6f\n",$7}' >> Gcorr

fi
done


index="name\t\t\t\tE0\t\tEzpe\t\tHcorr\t\tGcorr"
paste name E0 Ezpe Hcorr Gcorr | sed -r "1i$index" > energy.out
rm name E0 Ezpe Hcorr Gcorr

