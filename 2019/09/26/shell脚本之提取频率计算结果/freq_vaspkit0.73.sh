#!/bin/sh

touch run.sh
chmod +x run.sh

for dir in `find ./ -type d | grep 'freq'`; do
if [ ! -f "$dir/OUTCAR" ];then
continue
else
cat > run.sh <<EOR
#!/bin/sh
cd $dir
vaspkit <<!
5
501
$1
!
EOR

echo $dir  
./run.sh | tee out.log

echo $dir | sed -r 's/freq//' | sed 's/^..//' | awk '{printf "%-30s\n",$1}' >> name
grep entropy $dir/OUTCAR | tail -1 | awk '{printf "%.6f\n",$7}' >> E0
grep "Zero-point energy E_ZPE" out.log | awk '{printf "%.6f\n",$5}' >> Ezpe
grep "Enthalpy contribution E_H" out.log | awk '{printf "%.6f\n",$5}' >> H
grep "Entropy contribution T*S" out.log | awk '{printf "%.6f\n",$5}' >> TS
grep "Thermal correction to G(T)" out.log | awk '{printf "%.6f\n",$6}' >> G

fi
done


index="name\t\t\t\tE0\t\tEzpe\t\tH(T)\t\tG(T)"
paste name E0 Ezpe H TS G | sed -r "1i$index" > energy.out
rm name E0 Ezpe H TS G out.log run.sh

