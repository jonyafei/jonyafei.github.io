#!/bin/sh

touch run.sh
chmod +x run.sh

for dir0 in `find ./ -type d | grep 'opt'`; do
if [ ! -f "$dir0/OUTCAR" ];then
continue
else
#dir=`echo $dir0 | sed -r 's/.opt//'`
dir=${dir0%/opt}
#echo $dir | sed 's/^..//' | awk '{printf "%-30s\n",$1}' >> name
echo ${dir#*/} | awk '{printf "%-30s\n",$1}' >> name
grep entropy $dir/opt/OUTCAR | tail -1 | awk '{printf "%.6f\n",$7}' >> E0

if [ -f "$dir/freq/OUTCAR" ];then
cat > run.sh <<EOR
#!/bin/sh
cd $dir/freq
vaspkit <<!
5
501
$1
!
EOR

echo $dir 
./run.sh | tee out.log

grep "Zero-point energy E_ZPE" out.log | awk '{printf "%.6f\n",$7}' >> Ezpe
grep "Thermal correction to H(T)" out.log | awk '{printf "%.6f\n",$7}' >> H
grep "Thermal correction to G(T)" out.log | awk '{printf "%.6f\n",$7}' >> G

else
echo "None" | awk '{printf "%-8s\n",$1}' >> Ezpe
echo "None" | awk '{printf "%-8s\n",$1}' >> H
echo "None" | awk '{printf "%-8s\n",$1}' >> G

fi
fi
done


index="name\t\t\t\tE0\t\tEzpe\t\tH(T)\t\tG(T)"
paste name E0 Ezpe H G | sort | sed -r "1i$index" > energy.out
rm name E0 Ezpe H G out.log run.sh

