function [g Dvec up nd] =  vertSec(row,data)

xr = smooth(row,0.01,'loess');
   g=f_cumul(xr);
   
   Dvec = movingslope(g,2,1,data.dt);
   Dvec = Dvec - 0.05*max(Dvec); 
%       Dvec = Dvec - mean(Dvec); 
   Dvec(Dvec < 0 )= 0;
   [up nd] = araivEst(Dvec);
