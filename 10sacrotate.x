#!/bin/sh

echo "########    renaming the original sac files andd rotating the horizontal components    ##########"
echo "#################################################################################################"

e=` awk 'END{ print NR}' evfullist2.dat`

for ev in ` awk '{ print $1}' evfullist2.dat`; do
echo $ev

cd $ev
cd Icorrf/Vel

ls * | awk -F"." '{print $2,$3}' | sort | uniq -c > stn_3c.id

## rename the SAC filename to something like NET.STN.BH[ZNE].SAC
## find the station name which has three seismograms (three components)
#ls *.SAC | awk -F"." '{print $7,$8}' | sort | uniq -c | awk '{if ($1 == 3) print $2,$3}' | sort | uniq > stn_3c.id
## rename the SAC filename to NET.STN.BH[ZNE].SAC
#ls `awk '{print "*"$1"."$2".*SAC"}' stn_3c.id` | awk -F"." '{print "mv "$0,$7"."$8"."$10"."$12}' | sh
## move the rest into backup for now
#mv *.SAC backup

ls * > saclist5.dat
awk -F"." '{print $2}' saclist5.dat | sort -u | uniq -u > stnam5.dat
awk -F"." '{print $4}' saclist5.dat | sort -u | uniq -u > ty.dat

cp stnam5.dat stnam5.2.dat
k=` awk 'END{ print NR}' stnam5.2.dat`

echo -------------------------

for stn in ` awk '{ print $2}' stnam5.2.dat`; do

p=` awk '{print $1}' ty.dat`
echo $p
echo "r CN.$stn.??N.* " > 5macro.mac
echo "ch CMPAZ 0 CMPINC 90 " >> 5macro.mac
echo "wh " >> 5macro.mac
echo "r CN.$stn.??E.* " >> 5macro.mac
echo "ch CMPAZ 90 CMPINC 90 " >> 5macro.mac
echo "wh " >> 5macro.mac
echo "r CN.$stn.??N.* CN.$stn.??E.* " >> 5macro.mac
echo "rotate to GCP " >> 5macro.mac
echo "w CN.$stn.HHR.D.$p *.CN.$stn.HHT.D.$p " >> 5macro.mac
echo "q " >> 5macro.mac

sac 5macro.mac

tail -"$k" stnam5.2.dat > 3.2temp
mv 3.2temp stnam5.2.dat

done

# cut the two horizontal component data so that they will have the same length, which are needed for rotation in SAC.
# list the end of the two horizontal component data

#saclst e f `awk '{print $2"."$3".?H[EN].$p"}' stn_3c.id` 
#saclst e f `awk '{print $2"."$3".?H[EN].$p"}' stn_3c.id` |\
#paste - - | awk '{if ($2> $4) print $4;else print $2}' > BH_EN_end_time.dat
#saclst b f `awk '{print $2"."$3".?H[EN].$p"}' stn_3c.id` |\
#paste - - | awk '{if ($2< $4) print $4;else print $2}' > BH_EN_start_time.dat
# find out the minimum of the end time and save it
#paste -d" " stn_3c.id BH_EN_start_time.dat BH_EN_end_time.dat |\
#awk '{print "cut "$3,$4"\nr "$2"."$3".?H[EN].$p\nrmean\nrtrend\nrotate to GCP\nw "$2"."$3".#HR.$p "$2"."$3".#HT.$p"} END {print "q"}' | sac
# cut the two horizontal component data, remove mean and trend,
# rotate them to great circle path, save the output
# awk '{print "r "$1"."$2".*.$p\nrmean\nrtrend\nw over"} END {print "q"}' stn_3c.id | sac
# remove mean and trend for the original three component data.


#rm saclist5.dat
#rm -f stnam5.dat
#rm -f stnam4.dat
#rm -f stn_3c.id
#rm -f stnam5.2.dat
#rm -f stnam5.2.dat
#rm -f BH_EN_end_time.dat
#rm -f BH_EN_start_time.dat
#rm -f -r 5macro.mac
cd ../

cd ../../

tail -"$e" evfullist2.dat > 5temp
mv 5temp evfullist2.dat

done
rm -f 5temp

