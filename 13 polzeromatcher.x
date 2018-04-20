#!/bin/sh

# this script was written to find mach pole zero file of each sac files 


k=` awk 'END{ print NR}' stnam4.dat`

for stn in ` awk '{ print $1}' stnam4.dat`; do
arr=(` ls *.SAC | grep "$stn"`) 
n=${#arr[*]}

st=` echo ${arr[0]} | awk 'BEGIN{FS="."}{ print $8}'`
cm=` echo ${arr[0]} | awk 'BEGIN{FS="."}{ print $10}'`

st1=` echo ${arr[1]} | awk 'BEGIN{FS="."}{ print $8}'`
cm1=` echo ${arr[1]} | awk 'BEGIN{FS="."}{ print $10}'`

st2=` echo ${arr[2]} | awk 'BEGIN{FS="."}{ print $8}'`
cm2=` echo ${arr[2]} | awk 'BEGIN{FS="."}{ print $10}'`

echo ${arr[0]}

echo "r ${arr[0]}" > 5macro.m
echo "lh STLA STLO STEL" >> 5macro.m
echo "q" >> 5macro.m
sac 5macro.m > lh.s

stla=`sed -n '9p' lh.s | awk -F "=" '{print $2}'`
stlo=`sed -n '10p' lh.s | awk -F "=" '{print $2}'`
stel=`sed -n '11p' lh.s | awk -F "=" '{print $2}'`

stla1=`echo ${stla} | sed 's/[eE]+*/*10^/g' | bc -l`
stlo1=`echo ${stlo} | sed 's/[eE]+*/*10^/g' | bc -l`
stel1=`echo ${stel} | sed 's/[eE]+*/*10^/g' | bc -l`

echo $stla1
echo $stlo1
echo $stel1

for i in  $(ls SAC_PZs_CN_"$st"_"$cm"__*); do
echo $i

stla2=`grep LATITUDE $i | awk -F ":" '{print $2}'`
stlo2=`grep LONGITUDE $i | awk -F ":" '{print $2}'`
stel2=`grep ELEVATION $i | awk -F ":" '{print $2}'`

echo $stla2
echo $stlo2
echo $stel2


if [[ $(echo "if ($stla1 == $stla2) 1 else 0" | bc) -eq 1 ]];
 then

np=PZ.$st.$cm
cp $i $np
cp "SAC_PZs_CN_"$st"_"$cm1"__"* PZ.$st.$cm1
cp "SAC_PZs_CN_"$st"_"$cm2"__"* PZ.$st.$cm2
 
fi 

done
let k=$(($k-1))

tail -"$k" stnam4.dat > 4.2temp
mv 4.2temp stnam4.dat


done 
