#!/bin/bash

echo "###################################  7ampmaker.x  ################################################"
echo "########################        amplitute calculation       ######################################"
echo "##################################################################################################"

e=` awk 'END{ print NR}' evfullist2.dat`

for ev in ` awk '{ print $1}' evfullist2.dat`; do
echo -------------------------
echo " Number oof Events: "$e
echo -------------------------
echo "Event Name:"$ev
echo -------------------------
cd $ev
cd Icorrf
cd Vel

rm -f $ev.vel.out
ls CN.*.???.vel | awk -F"." '{print $2}' | sort | uniq -c > stn.dat
cp stn.dat ../.


for stn in $(ls *.??Z.vel); do
echo $stn

echo "r $stn" > 3.1.macro.mac
echo "lh OMARKER DIST" >> 3.1.macro.mac
echo "q" >> 3.1.macro.mac
sac 3.1.macro.mac > o.dis.lh

origin=`sed -n '9p' o.dis.lh | awk -F "=" '{print $2}'`
dist=`sed -n '10p' o.dis.lh | awk -F "=" '{print $2}'`

o=`echo ${origin} | sed 's/[eE]+*/*10^/g' | bc -l`
d=`echo ${dist} | sed 's/[eE]+*/*10^/g' | bc -l`
echo "origin time: "$o 
echo "distance: "$d

t1=$(echo "$d/5.8"+"$o" | bc -l)
t2=$(echo "$d/5"+"$o" | bc -l)

t3=$(echo "$d/3.65"+"$o" | bc -l)
t4=$(echo "$d/2.6"+"$o" | bc -l)

echo $t1
echo $t2 

echo $t3
echo $t4

st2=` echo ${stn} | awk 'BEGIN{FS="."}{ print $2}'`
cm2=` echo ${stn} | awk 'BEGIN{FS="."}{ print $3}'`
echo $st2
echo $cm2

echo "r ${stn}" > macrocut9.mac
echo "sync" >> macrocut9.mac
echo "MTW $t1 $t2 " >> macrocut9.mac
echo "RMS to user1 " >> macrocut9.mac
echo "MTW $t3 $t4 " >> macrocut9.mac
echo "RMS to user2 " >> macrocut9.mac
echo "lh user1 user2 " >> macrocut9.mac
echo "q" >> macrocut9.mac
sac  macrocut9.mac > u1.u2.lh

us21=`sed -n '9p' u1.u2.lh | awk -F "=" '{print $2}'`
us22=`sed -n '10p' u1.u2.lh | awk -F "=" '{print $2}'`
echo "us21"$us21
echo $us22

u21=`echo ${us21} | sed 's/[eE]+*/*10^/g' | bc -l`
u22=`echo ${us22} | sed 's/[eE]+*/*10^/g' | bc -l`

echo "u21"$u21 
echo $u22

snr2=$(echo "$u22/$u21" | bc -l)
echo "snr "$snr2

echo "r ${stn}" > macrocut10.mac
echo "cut $t3 $t4" >> macrocut10.mac
echo "r ${stn}" >> macrocut10.mac
echo "w $st2"."$cm2".lg.vel"" >> macrocut10.mac
echo "q" >> macrocut10.mac
sac macrocut10.mac


#arr=(` echo 0.3 1.3 0.5 1.5 0.9 1.9 1.5 2.5 2.3 3.3 3.3 4.3 4.5 5.5 5.9 6.9 7.5 8.5 9.3 10.3 11.3 12.3 13.5 14.5` )
arr=(` echo 0.5 1.0 1.0 2.0 2.0 4.0 4.0 8.0 8.0 16.0` )

n=${#arr[*]}

cp ../*.SAC ./

#for i in $(ls *.lg.vel); do
echo "r $st2"."$cm2".lg.vel"" > seis2.macro
echo "lh KNETWK KSTNM KCMPNM EVLO EVLA EVDP STLO STLA MAG DIST DELTA" >> seis2.macro
echo "q" >> seis2.macro

echo $ev | awk 'BEGIN{FS="/"}{printf( "%s ", $1); }' >> $ev.vel.out
sac seis2.macro | grep "KNETWK" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "KSTNM" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "KCMPNM" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "EVLO" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "EVLA" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "EVDP" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "STLO" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "STLA" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "MAG" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "DIST" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
sac seis2.macro | grep "DELTA" | awk '{printf( "%s ", $3); }' >> $ev.vel.out
#ls *.SAC | awk -F"." -v s="$snr2" '{printf( "%s ", $1" " $2" " $3" " $4" " $5" " $6" 0 " s);}' >> $ev.vel.out
ls *.SAC | awk -F"." '{printf( "%s ", $1" " $2" " $3" " $4" " $5" " $6" 0 1000 ");}' >>$ev.vel.out


for  ((j=0;j<=n-1;j+=2)); do
echo "r $st2"."$cm2".lg.vel"" > seis1.macro
#echo "MTW $t1 $t2 " >> seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[j]} ${arr[j+1]}" >> seis1.macro
echo "fft " >> seis1.macro
echo "keepam" >> seis1.macro
echo "rms to user3 " >> seis1.macro
echo "lh user3 " >> seis1.macro
echo "q" >> seis1.macro
sac seis1.macro | grep "user3" | awk '{if (NF == 3) printf( "%s ", $3);}' >> $ev.vel.out

done  

echo " " | awk '{printf(" \n");}' >> $ev.vel.out


done 
mv $ev.vel.out ../.

rm seis1.macro
rm seis2.macro
rm -f saclist.dat
rm -f stnam.d
rm -f 3stnam.dat
rm -f 3.1.macro.mac
rm -f o.dis.lh
rm -f 3.2.macro.mac
rm -f u1.u2.lh
rm -f 3.2temp
rm -f stn.dat
rm -f *.SAC
rm -f -r macrocut9.mac
rm -f -r macrocut10.mac

cd ..

cd Dsp

for stn in $(ls *.??Z.dsp); do
echo $stn


echo "r $stn" > 3.1.macro.mac
echo "lh OMARKER DIST" >> 3.1.macro.mac
echo "q" >> 3.1.macro.mac
sac 3.1.macro.mac > o.dis.lh

origin=`sed -n '9p' o.dis.lh | awk -F "=" '{print $2}'`
dist=`sed -n '10p' o.dis.lh | awk -F "=" '{print $2}'`

o=`echo ${origin} | sed 's/[eE]+*/*10^/g' | bc -l`
d=`echo ${dist} | sed 's/[eE]+*/*10^/g' | bc -l`
echo "origin time: "$o 
echo "distance: "$d

t1=$(echo "$d/5.8"+"$o" | bc -l)
t2=$(echo "$d/5"+"$o" | bc -l)

t3=$(echo "$d/3.65"+"$o" | bc -l)
t4=$(echo "$d/2.6"+"$o" | bc -l)

echo $t1
echo $t2 

echo $t3
echo $t4

st2=` echo ${stn} | awk 'BEGIN{FS="."}{ print $2}'`
cm2=` echo ${stn} | awk 'BEGIN{FS="."}{ print $3}'`

echo $cm2

echo "r $stn" > macrocut9.mac
echo "sync" >> macrocut9.mac
echo "MTW $t1 $t2 " >> macrocut9.mac
echo "RMS to user1 " >> macrocut9.mac
echo "MTW $t3 $t4 " >> macrocut9.mac
echo "RMS to user2 " >> macrocut9.mac
echo "lh user1 user2 " >> macrocut9.mac
echo "q" >> macrocut9.mac
sac  macrocut9.mac > u1.u2.lh

us21=`sed -n '9p' u1.u2.lh | awk -F "=" '{print $2}'`
us22=`sed -n '10p' u1.u2.lh | awk -F "=" '{print $2}'`

u21=`echo ${us21} | sed 's/[eE]+*/*10^/g' | bc -l`
u22=`echo ${us22} | sed 's/[eE]+*/*10^/g' | bc -l`

snr2=$(echo "$u22/$u21" | bc -l)

echo "r $stn" > macrocut10.mac
echo "cut $t3 $t4" >> macrocut10.mac
echo "r $stn" >> macrocut10.mac
echo "w $st2"."$cm2".lg.dsp"" >> macrocut10.mac
echo "q" >> macrocut10.mac
sac macrocut10.mac

#arr=(` echo 0.3 1.3 0.5 1.5 0.9 1.9 1.5 2.5 2.3 3.3 3.3 4.3 4.5 5.5 5.9 6.9 7.5 8.5 9.3 10.3 11.3 12.3 13.5 14.5` )
arr=(` echo 0.5 1.0 1.0 2.0 2.0 4.0 4.0 8.0 8.0 16.0` )

n=${#arr[*]}

cp ../*.SAC ./

echo "r $st2"."$cm2".lg.dsp "" > seis2.macro
echo "lh KNETWK KSTNM KCMPNM EVLO EVLA EVDP STLO STLA MAG DIST DELTA" >> seis2.macro
echo "q" >> seis2.macro

echo $ev | awk 'BEGIN{FS="/"}{printf( "%s ", $1); }' >> $ev.dsp.out
sac seis2.macro | grep "KNETWK" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "KSTNM" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "KCMPNM" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "EVLO" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "EVLA" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "EVDP" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "STLO" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "STLA" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "MAG" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "DIST" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "DELTA" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
ls *.SAC | awk -F"." '{printf( "%s ", $1" " $2" " $3" " $4" " $5" " $6" 0 1000 ");}' >>$ev.dsp.out

for  ((j=0;j<=n-1;j+=2)); do
echo "r $st2"."$cm2".lg.dsp"" > seis1.macro
#echo "MTW $t1 $t2 " >> seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[j]} ${arr[j+1]}" >> seis1.macro
echo "fft " >> seis1.macro
echo "keepam" >> seis1.macro
echo "rms to user4" >> seis1.macro
echo "lh user4 " >> seis1.macro
echo "q" >> seis1.macro
sac seis1.macro | grep "user4" | awk '{if (NF == 3) printf( "%s ", $3);}' >> $ev.dsp.out
done 

echo " " | awk '{printf(" \n");}' >> $ev.dsp.out

done 
mv $ev.dsp.out ../.

rm seis1.macro
rm seis2.macro
rm -f saclist.dat
rm -f stnam.d
rm -f 3stnam.dat
rm -f 3.1.macro.mac
rm -f o.dis.lh
rm -f 3.2.macro.mac
rm -f u1.u2.lh
rm -f 3.2temp
rm -f stn.dat
rm -f *.SAC
rm -f -r macrocut9.mac
rm -f -r macrocut10.mac

cd ..
cd ../../

done



