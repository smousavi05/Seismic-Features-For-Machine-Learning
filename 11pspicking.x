
#!/bin/sh

echo "########################################  11pspicking.x  ########################################"
echo "#################       picking P and S arrivals  and making quality.dat     ####################"
echo "#################################################################################################"

e=` awk 'END{ print NR}' evfullist2.dat`

for ev in ` awk '{ print $1}' evfullist2.dat`; do
echo -------------------------
echo " Number oof Events: "$e
echo -------------------------
echo "Event Name:"$ev
echo -------------------------
cd $ev
cd CNDC
                                   
awk -F"." '{print $8}' saclist.dat | sort -u | uniq -u > stnam.d

echo "Record Quality: " > quality.dat
echo "A - Good" >> quality.dat
echo "B - Bad, problem with component" >> quality.dat
echo "? - Accepted, but might have problem" >> quality.dat
echo "C - Clipped" >> quality.dat
echo "D - directivity" >> quality.dat
echo "N - radiation effect" >> quality.dat
echo "M - missing data in waveform" >> quality.dat
echo "O - Orientation missing" >> quality.dat
echo "S - spikes in data" >> quality.dat
echo "X - no data - in noise" >> quality.dat
echo "net stn cmp  comment" >> quality.dat

cp stnam.d 3stnam.dat
k=` awk 'END{ print NR}' 3stnam.dat`
echo "Number of Stations:"$k
echo -------------------------
for stn in ` awk '{ print $1}' 3stnam.dat`; do
ls *.SAC | grep "$stn" > 3temp
arr=(` ls *.SAC | grep "$stn"`) 
n=${#arr[*]}
if [ "$n" == "3" ]; then
echo $n
echo ${arr[*]}
echo "r ${arr[0]} ${arr[1]} ${arr[2]}" > 3macro.m
echo "qdp off " >> 3macro.m
echo "rmean" >> 3macro.m
echo "sync" >> 3macro.m
echo "hp bu co 0.5 n 4 p 2 " >> 3macro.m
echo "p1" >> 3macro.m
sac 3macro.m 
echo ---------------------------------------
echo "A - Good" 
echo "B - Bad, problem with component" 
echo "? - Accepted, but might have problem" 
echo "C - Clipped" 
echo "D - directivity" 
echo "N - radiation effect" 
echo "M - missing data in waveform" 
echo "O - Orientation missing" 
echo "S - spikes in data" 
echo "X - no data - in noise" 
echo ----------------------------------------
echo what you think A B ? C D N M O S X
read com
                         
awk -F"." -v c="$com" '{print $7".."$8"."$10"  "c}' 3temp >> quality.dat

if [ "$com" == "A" ]; then 
echo "r ${arr[0]} ${arr[1]} ${arr[2]}" > 3macro.m
echo "qdp off " >> 3macro.m
echo "rmean" >> 3macro.m
echo "sync" >> 3macro.m
echo "bp bu p 2 n 4 c 0.06 0.09" >> 3macro.m
echo "p1" >> 3macro.m
echo "lh DIST AMARKER T0MARKER GCARC EVDP" >> 3macro.m
echo  -------------------------------------------------
echo " Now Pickk the P & S ( hint: xlim  ppk m   wh ) "
echo  -------------------------------------------------
sac 3macro.m 

elif [ "$com" == "?" ]; then 
echo "r ${arr[0]} ${arr[1]} ${arr[2]}" > 3macro.m
echo "qdp off " >> 3macro.m
echo "rmean" >> 3macro.m
echo "sync" >> 3macro.m
echo "bp bu p 2 n 4 c 0.06 0.09" >> 3macro.m
echo "p1" >> 3macro.m
echo "lh DIST AMARKER T0MARKER GCARC EVDP" >> 3macro.m
echo  -------------------------------------------------
echo " Now Pickk the P & S ( hint: xlim  ppk m   wh ) "
echo  -------------------------------------------------
sac 3macro.m 
fi 

let k=$(($k-1))
echo " Remaining Stations:"$k

tail -"$k" 3stnam.dat > 3.2temp
mv 3.2temp 3stnam.dat
read -p next
rm -f 3macro.m
fi 
 
done 
mv quality.dat ..

cd ../../
let e=$(($e-1))

tail -"$e" evfullist2.dat > 3.3temp
mv 3.3temp evfullist2.dat
read -p next
rm -f 3.3temp

done
