#!/bin/sh

echo "#########################################  fft and noise    #####################################"
echo "#################################################################################################"

e=` awk 'END{ print NR}' evfullist2.dat`

for ev in ` awk '{ print $1}' evfullist2.dat`; do
echo -------------------------
echo " Number oof Events: "$e
echo -------------------------
echo "Event Name:"$ev
echo -------------------------
cd $ev
cd Icorrf

echo "cut 0 A " > noise.mac
echo "r Vel/* " >> noise.mac
echo "cut off " >> noise.mac
echo "w change vel noise.vel " >> noise.mac
echo "q " >> noise.mac

#sac noise.mac

echo "r Vel/* " > fft.mac
echo "fft " >> fft.mac
echo "wsp am " >> fft.mac
echo "q " >> fft.mac

#sac fft.mac

cd Vel

ls *.vel | awk -F"." '{print $2,$3}' | sort | uniq -c > stn.dat
k=` awk 'END{ print NR}' stn.dat`
echo "Number of Stations:"$k
echo -------------------------
for sc in ` awk '{ print $2"."$3}' stn.dat`; do

echo "qdp 70000 " > 6macro.mac
echo "loglog" >> 6macro.mac
echo "xlim .01 50 " >> 6macro.mac
echo "color 1 i l 7 1 " >> 6macro.mac
echo "r *.$sc.*.am " >> 6macro.mac
echo "p2 " >> 6macro.mac
echo "ppk " >> 6macro.mac

sac 6macro.mac


tail -"$k" stn.dat > temp
mv temp stn.dat
read -p next

done 
cd ..

rm -f 3.2temp
#rm -f -r noise.mac
#rm -f -r ftt.mac

cd ../../

tail -"$e" evfullist2.dat > 3.2temp
mv 3.2temp evfullist2.dat

done 

