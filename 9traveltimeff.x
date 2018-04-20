#!/bin/bash

echo "################################  9traveltimeff.x  ##########################################" 
echo "########################        Travel tim plot       #######################################"
echo "#############################################################################################"

echo "which azimuth ?"
read az 
echo "which velocity ?"
read vel

for ev in ` awk '{ print $1}' evfullist2.dat`; do
echo -------------------------
echo "Event Name:"$ev
echo -------------------------
cd $ev
cd Icorrf
cd Dsp

rm -f prs.macro

rm -f $az.list
rm -f azdistlist.txt 

azb=$(echo "$az-15" | bc -l)
aze=$(echo "$az+15" | bc -l)

echo "lower limit of azimouth: "$azb
echo "upper limit of azimouth: "$aze

for i in ` ls *.dsp `; do

echo $i

echo "r $i" > rddisaz.macro
echo "lh KSTNM AZ BAZ  " >> rddisaz.macro
echo "q" >> rddisaz.macro

sac rddisaz.macro | grep "KSTNM" | awk '{printf( "%s", $3" "); }' >> azdistlist.txt 
sac rddisaz.macro | grep "AZ" | awk '{printf( "%s", $3" "); }' >> azdistlist.txt
sac rddisaz.macro | grep "BAZ" | awk '{printf( "%s", $3" "); }' >> azdistlist.txt
echo " " | awk '{printf(" \n");}' >> azdistlist.txt

done 

for st in ` awk '{ print $1}'  azdistlist.txt`; do
echo $st
grep "$st"  azdistlist.txt > temp11

a=` cat temp11 | awk -F" " '{printf( "%s ",$2);}'`
echo $a

azz=`echo ${a} | sed 's/[eE]+*/*10^/g' | bc -l`
echo $azz

if [[ $(echo "if ($azz >= $azb) 1 else 0" | bc) -eq 1 && $(echo "if ($azz <= $aze) 1 else 0" | bc) -eq 1 ]]; then

#echo "$st" >> $az.list

#cat $az.list | awk '{printf ("%s ", " r more CN."$st".??Z.dsp");}' >> prs.macro

echo " r more "*."$st".*.dsp" >> prs.macro

echo " " | awk '{printf(" \n");}' >> prs.macro

fi 

done 

#if [[ -a $az.list ]]; then


echo "INTERPOLATE DELTA 0.01 " >> prs.macro
echo "qdp off " >> prs.macro

echo "sss" >> prs.macro

echo "vm 1 refractedwave vapp $vel" >> prs.macro
echo "timewindow 0 600" >> prs.macro
echo "title on "Event: $ev "" >> prs.macro  
echo "bd sgf" >> prs.macro 
echo "prs O r" >> prs.macro
echo "ed sgf" >> prs.macro 
echo "sgftops f001.sgf f001.ps" >> prs.macro
echo " q" >> prs.macro

sac prs.macro
mv f001.ps $az.$vel.$ev.ps
gv $az.$vel.$ev.ps

#mv $az.$vel.$ev.ps ../.

#else 
#echo  "####################################################################"
#echo "# No stations within this azimouth range in $ev # " 
#echo  "####################################################################"

#fi

rm -f f001.sgf
#rm -f azdistlist.txt
rm -f temp11
#rm -f prs.macro
rm -f rddisaz.macro

cd ../../../
done 
