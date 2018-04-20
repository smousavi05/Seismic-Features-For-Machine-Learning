#!/bin/bash

echo "##########################################################################################"
echo "########################        amp maker                      ######################################"
echo "##################################################################################################"

for ev in ` awk '{ print $1}' f.dat`; do
echo $ev
cd $ev
rm -f *.out

########
for stn in ` ls *,dsp`;do

# reading the distance and origin time and calculating Lg and p coda's window and calculating RMS and SNR

echo "r $stn" > 3.1.macro.mac
echo "lh OMARKER DIST" >> 3.1.macro.mac
echo "q" >> 3.1.macro.mac
sac 3.1.macro.mac > o.dis.lh

origin=`sed -n '9p' o.dis.lh | awk -F "=" '{print $2}'`
dist=`sed -n '10p' o.dis.lh | awk -F "=" '{print $2}'`

o=`echo ${origin} | sed 's/[eE]+*/*10^/g' | bc -l`
d=`echo ${dist} | sed 's/[eE]+*/*10^/g' | bc -l`

######################################
t1=$(echo "$d/5.8"+"$o" | bc -l)   # Pn
t2=$(echo "$d/4.8"+"$o" | bc -l)

#t3=$(echo "$d/3.7"+"$o" | bc -l) # Lg
#t4=$(echo "$d/3.0"+"$o" | bc -l)

t3=$(echo "$d/4.8"+"$o" | bc -l) # Sn
t4=$(echo "$d/3.6"+"$o" | bc -l)

t5=$(echo "$d/3.5"+"$o" | bc -l)
#######################################
echo "r $stn" > 3.2.macro.mac
echo "sync" >> 3.2.macro.mac
echo "MTW $t1 $t2 " >> 3.2.macro.mac
echo "RMS to user1 " >> 3.2.macro.mac
echo "MTW $t3 $t4 " >> 3.2.macro.mac   
echo "RMS to user2 " >> 3.2.macro.mac
echo "lh user1 user2 " >> 3.2.macro.mac
echo "q" >> 3.2.macro.mac
sac 3.2.macro.mac > u1.u2.lh

st=` echo $stn | awk 'BEGIN{FS="."}{ print $3}'`
us01=` sed -n '9p' u1.u2.lh | awk -F "=" '{print $2}'`
us02=` sed -n '10p' u1.u2.lh | awk -F "=" '{print $2}'`

u01=`echo ${us01} | sed 's/[eE]+*/*10^/g' | bc -l`
u02=`echo ${us02} | sed 's/[eE]+*/*10^/g' | bc -l`
snr0=$(echo "$u02/$u01" | bc -l)
SNR=`echo ${snr0} | sed 's/[eE]+*/*10^/g' | bc -l`
echo $SNR

############################################

st=` echo $stn | awk 'BEGIN{FS="-"}{ print $3}'`
y=` echo $stn | awk 'BEGIN{FS="-"}{ print $1}'`
t=` echo $stn | awk 'BEGIN{FS="-"}{ print $2}'`

echo "***************"
echo "station "$st
echo "origin time: "$o 
echo "distance: "$d
echo "tlgb" $t3
echo "tlge" $t4
echo "t5" $t5

echo "***************"

# taking sac header info and putting into output file 
echo "r $stn" > seis2.macro
echo "lh KNETWK KSTNM KCMPNM EVLO EVLA EVDP STLO STLA MAG DIST DELTA KZDATE KZTIME" >> seis2.macro
echo "q" >> seis2.macro

echo $stn | awk 'BEGIN{FS="-"}{printf( "%s ", $1"."$2); }' >> $ev.dsp.out
sac seis2.macro | grep "EVLA" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "EVLO" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "KSTNM" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "STLA" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "STLO" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
sac seis2.macro | grep "DIST" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out
#echo $t5 | awk -v s="$t5" '{printf( "%s ",s);}' >> $ev.dsp.out
sac seis2.macro | grep "MAG" | awk '{printf( "%s ", $3); }' >> $ev.dsp.out

echo $stn | awk 'BEGIN{FS="-"}{ printf( "%s ", $4); }'  >> $ev.dsp.out
echo $stn | awk 'BEGIN{FS="-"}{ printf( "%s ", $5); }' | awk 'BEGIN{FS=","}{ printf( "%s", $1); }' >> $ev.dsp.out

echo $SNR | awk '{ printf( "%s "," "$1); }'  >> $ev.dsp.out

# filtering the seismoram and writing down the spectogram
echo "r $stn" > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "MTW $t3 $t4 " >>  seis1.macro
echo "fft " >> seis1.macro
#echo "sync " >> seis1.macro
echo "WRITESP am $stn".fft"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

# calculating RMS of spectral amplitude
echo "r $stn.fft.am" > seis3.macro  
echo "smooth " >> seis3.macro
echo "RMS to user3 " >>  seis3.macro
echo "lh user3 " >>  seis3.macro
echo "q " >>  seis3.macro
sac seis3.macro > u2.lh
sed -n '9p' u2.lh | awk -F "=" '{printf("%s",$2);}'>> $ev.dsp.out

# calculating the maximum amplitude of hilbert transform envelop
#echo "read $st2"."$cm2".fft.dsp.am"" > seis3.macro  
#echo "smooth " >> seis3.macro
#echo "envelope">> seis3.macro
#echo "lh DEPMAX " >> seis3.macro
#echo "q" >> seis3.macro
#sac seis3.macro | grep "DEPMAX" | awk '{if (NF == 3) printf( "%s ", $3);}'>> $ev.dsp.out

echo " " | awk '{printf(" \n");}' >> $ev.dsp.out

done 

rm seis1.macro
rm seis2.macro
rm seis3.macro
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

done 


