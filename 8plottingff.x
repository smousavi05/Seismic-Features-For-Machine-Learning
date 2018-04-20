#!/bin/sh

echo "#####################################  8plottingff.x  ###########################################"
echo "#################################         plotting         ######################################"
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

echo ">" > path.gmt
echo "" > station.gmt

echo 
for stn in ` awk '{ print $2}' stn.dat`; do

grep "$stn" $ev.vel.out > temp
cat temp | awk -F" " '{printf( "%s ",$5" "$6);}' >> path.gmt
echo " " | awk '{printf(" \n");}' >> path.gmt
cat temp | awk -F" " '{printf( "%s ",$8" "$9);}' >> path.gmt
echo " " | awk '{printf(" \n");}' >> path.gmt
echo ">" >> path.gmt

cat temp | awk -F" " '{printf( "%s ",$8" "$9);}' >>  station.gmt
echo " " | awk '{printf(" \n");}' >> station.gmt

done

REGION=-90/-45/35/60
size=20c
fromXaxis=4
fromYaxis=3
title=${ev}

##########   set defaults
WHITE=255
LTGRAY=192
VLTGRAY=225
EXTGRAY=250
GRAY=128
BLACK=0
RED=250/0/0
DKRED=196/50/50
BLUE=0/0/255
GREEN=0/255/0
YELLOW=255/255/50
ORANGE=255/192/50
PURPLE=255/50/255
CYAN=50/255/255
LTBLUE=192/192/250
VLTBLUE=225/250/250
LTRED=250/225/225
PINK=255/225/255
BROWN=160/64/32

gmtset MEASURE_UNIT cm

########  basemap
psbasemap -R${REGION} -JM${size} -Ba2f2g1/a2f2g1:.${title}: -X${fromXaxis} -Y${fromYaxis} -K > $0.ps

########  converting the img file to the grd file
img2grd ../../../topo_8.2.img -R$REGION -Gtopo2.grd -T1 

########  making the CPT
makecpt -Cglobe -T-12000/12000/500 -Z > topo.cpt

########   grdgradient grd to gradient
grdgradient ../../../topo2.grd -Gtopo.grad -A5 -Nt0.9

##########  ploting the grd file on map
grdimage ../../../topo2.grd -Ctopo.cpt -Itopo.grad -JM${size} -R${REGION} -B1000g10000nSeW -Sb -K -O >> $0.ps

##########   cost line data
pscoast -R${REGION} -JM${size} -W1 -Df+ -A10000 -K -O -Na/3/0 >> $0.ps

##### DRAW LINES 
psxy path.gmt -R${REGION} -JM${size} -M -K -O -W1  >> $0.ps

###### plot stations
cat station.gmt | awk -F" " '{print("\n", $1" "$2);}' > lo.xy
psxy lo.xy -R${REGION} -JM${size} -St.2 -M -W1/0 -G$RED -K -O >> $0.ps

###### plot Earthquake 
awk '{print $5,$6}' temp | head -1000 > locations.xy
psxy locations.xy -R${REGION} -JM${size} -Sc0.2c -W1/0 -G$YELLOW -M -K -O >> $0.ps

########## Scale Bar
pscoast -R${REGION} -JM${size} -W1 -K -O -Lf-93/32.5/32.7/100 -V >> $0.ps 

#########    plotting legend
pstext -R -J -O -K << END >> $0.ps
-93.5 33.4 12 0 0 0 Stations
-93.5 33.2 12 0 0 0 Earthquake	
END

echo -89.0 36.45 | psxy -R -J -St.2 -W1/0 -G$RED -O -K >> $0.ps
echo -89.0 36.25 | psxy -R -J -Sc0.2c -W1/0 -G$YELLOW -O  >> $0.ps

#gv $0.ps &


rm topo2.grd 
rm topo.grad
rm topo.cpt
#rm station.gmt
rm lo.xy
rm locations.xy
#rm temp
#rm evlist.dat
#rm path.gmt

cd ../../
done