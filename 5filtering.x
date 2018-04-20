#!/bin/bash

echo "##################################################################################################"
echo "########################        filtrring      ######################################"
echo "##################################################################################################"

rm -f [12345]/*

e=` awk 'END{ print NR}' evfullist2.dat`

for ev in ` awk '{ print $1}' evfullist2.dat`; do
echo "************************************"
echo " Number oof Events: "$e
echo "************************************"
echo "Event Name:"$ev
echo "************************************"
cd $ev
rm -f *.out
cd Icorrf

#########
cd Dsp

########
arr=(` echo 0.9 1.1 1.1 1.4 1.4 1.85 1.85 2.525 2.525 3.5735 3.5735 5.05625 5.05625 7.334375 7.334375 10.7515625` )
#arr=(` echo 1 1.5 1.5 2.25 2.25 3.375 3.375 5.065 5.065 7.6 7.6 11.38 11.38 17.07 ` )
#arr=(` echo 0.5 1.0 1.0 1.75 1.75 2.875 2.875 4.565 4.565 7.1 7.1 10.88 10.88 16.57 ` )
#arr=(` echo 0.5 1.0 1.0 2.0 2.0 4.0 4.0 8.0 8.0 16.0` )
#arr=(` echo 1 1.5 1.5 2.5 2.5 4.5 4.5 8.5 8.5 16.5` )

echo ${arr[*]}
n=${#arr[*]}

########
for stn in $(ls *.dsp); do

st=` echo ${stn} | awk 'BEGIN{FS="."}{ print $2}'`
d=` echo ${ev} | awk 'BEGIN{FS="_"}{ print $2}'| awk 'BEGIN{FS="."}{ print $1}'`
t=` echo ${ev} | awk 'BEGIN{FS="_"}{ print $2}'| awk 'BEGIN{FS="."}{ print $2}'`

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[0]} ${arr[1]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/1/$d"-"$t"-"$st"-"${arr[0]}"-"${arr[1]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[2]} ${arr[3]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/2/$d"-"$t"-"$st"-"${arr[2]}"-"${arr[3]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[4]} ${arr[5]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/3/$d"-"$t"-"$st"-"${arr[4]}"-"${arr[5]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[6]} ${arr[7]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/4/$d"-"$t"-"$st"-"${arr[6]}"-"${arr[7]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[8]} ${arr[9]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/5/$d"-"$t"-"$st"-"${arr[8]}"-"${arr[9]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[10]} ${arr[11]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/6/$d"-"$t"-"$st"-"${arr[10]}"-"${arr[11]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[12]} ${arr[13]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/7/$d"-"$t"-"$st"-"${arr[12]}"-"${arr[13]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

echo " r $stn " > seis1.macro
echo "rmean " >> seis1.macro
echo "rtrend" >> seis1.macro
echo "taper" >> seis1.macro
echo "bp bu n 4 p 2 co ${arr[14]} ${arr[15]}" >> seis1.macro
echo "sync " >>  seis1.macro
echo "WRITE ../../../filtered/8/$d"-"$t"-"$st"-"${arr[14]}"-"${arr[15]}",dsp"" >> seis1.macro 
echo "q" >> seis1.macro
sac seis1.macro

done 

cd ..
cd ../../

done

