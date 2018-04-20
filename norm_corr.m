function [c_coef, c_zero, tshift]=norm_corr(s1,s2);
%**************************************************************************
%
%  Correlate 2 time series.  
%
%  s1 = sacfile data structure
%  s2 = sacfile data structure
%  
%  c_coef = maximum of the normalized cross correlation
%  c_zero = value of the correlation at zero time shift
%  tshift = time shift of the maximum normalized correlation
%
%**************************************************************************

%  Make sure the sampling interval is the same. If not, punt.
if s1(1).samprate ~= s2(1).samprate;
    c_coef=0.0;
    c_zero=0.0;
    tshift=0.0;
    fprintf('********Error from norm_corr.  Sampling rate inconsistent! Skipping**********');
else
    
%  Set length of arrays to the smaller of data and synthetic
npvd=s1(1).nsamps;
npvs=s2(1).nsamps;

np=min(npvd,npvs);
dat_s1(1:np)=s1(1).data(1:np);
dat_s2(1:np)=s2(1).data(1:np);

dt=1.0./s1(1).samprate;
tplot=linspace(0,(np-1).*dt,np);

t=linspace(-dt.*(np),dt.*np,2*np-1);
nzero=np;

cross_cor=xcorr(dat_s1,dat_s2,'coeff');
[ccor_max, indx]=max(cross_cor);

c_coef=ccor_max;
c_zero=cross_cor(nzero);
tshift=(indx - nzero).*dt;

PLOTFLAG=0;

if PLOTFLAG == 1;
    figure;
    subplot(3,1,1);
    plot(tplot,dat_s1,'-k');
    subplot(3,1,2);
    plot(tplot,dat_s2,'-k');
    subplot(3,1,3);
    t1=t(indx);
    plot(t,cross_cor,'-k',t1,c_coef,'ok');
end;

end;





