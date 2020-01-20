
 
 load('/deepvv2/74/z.mat'); masterD=z;clear z;dt=1/200; % 74 507
 load('/deepvv2/507/z.mat'); masterS=z;clear z;
 
 % Parameters for Calculate the wavelet transform -
opt.type = 'morlet';         % Mother wavelet type
opt.padtype = 'symmetric';   % padded via symmetrization
opt.rpadded = 1;
opt.nv = 16; 
data.dt=1/200;
dt=1/200;
 
no=[];dp=[];la=[];lo=[];
 
maxA=[];  % Average of maximum amplitude on thhree components
SpecCent = []; % Average of Spectral centroid on three component
Slope = [];    % Slope 

rmsA = [];            % RMS of frequency amplitude 
maxPFA = [];       % Maximum Power of Frequency Amplitude
maxFA = [];          % Dominent frequency
maxPF_FA = []; % Maximum power of frequency amplitude / dominent frequency
skwnss = [];
 
Energy = []; % energy 
 
otY = [];otX = [];otZ = [];
 
rec = []; azimuth = []; via = [];  DIP = []; dd_re = []; DIPp = [];
    
dr=[];ccD=[];ccS=[];xpD=[];xpS=[];xppD=[];xppS=[];

ccnAb2D=[]; ypicD=[]; ccnRel2D=[]; ccAb2D=[]; ccRel2D=[];
ccnAb2S=[]; ypicS=[]; ccnRel2S=[]; ccAb2S=[]; ccRel2S=[];

energy = [];
  
xppD=[]; xppS=[]; xpp2D=[]; xpp2S=[];
semSD=[];semSS=[];semDD=[];semDS=[]; 

semblanceD = [];
semblanceS = [];

envelopD = [];
envelopS = [];
duD = [];
duS = [];
d2uD = [];
d2uS = [];

% figure
% plot(zz,d2uS,'o');
% hold on
% plot(zz(1:218),d2uS(1:218),'s','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',10)
% plot(zz(219:end),d2uS(219:end),'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8)
% xlabel('Event Number','FontSize',12,'FontWeight','bold')
% ylabel('envelop','FontSize',12,'FontWeight','bold')
% grid on
% axis tight 
% legend('', 'Deep','Shallow','FontSize',12)
 

 for ii=1:1033
    disp(ii)
    
    v = sprintf('/%s/%s/%s','deepvv2', num2str(ii),'sac');
    [t, d, hdr] = fget_sac(sprintf(v));
    
%     %%% header's info 
%     dp=[dp, hdr.event.evdp];
%     la=[la, hdr.event.evla];
%     lo=[lo, hdr.event.evlo];
%     no=[no, ii];
%     
   load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'zd.mat'));
   load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'xd.mat'));
   load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'yd.mat'));
   
   %%% maximum amplitude
   llz = abs(zd);[llz idz]= max(llz);
   llx = abs(xd);[llx idx]= max(llx);
   lly = abs(yd);[lly idy]= max(lly);
   ll = (llz + llx + lly)/3;
   maxA=[maxA ll];
  
   %%% slope
   d1z=zd(200); d1x=xd(200); d1y=yd(200);
   sz = (llz-d1z)./(idz - 200); 
   sx = (llx-d1x)./(idx - 200); 
   sy = (lly-d1y)./(idy - 200); 
   
   sl = (sz + sx + sy)/3;
   Slope = [Slope sl];
   
   
   %%% Spectral Centroid 
    [CZ, CZM, CSD, CMax] = SpecCentroid(zd,200,528,50,30);
    [CX, CXM, CSD, CMax] = SpecCentroid(xd,200,528,50,30);
    [CY, CYM, CSD, CMax] = SpecCentroid(yd,200,528,50,30);
  
  
%    spCen = (CZ + CX + CY)/3;
   spCen = (CZM + CXM + CYM)/3;
   
   SpecCent = [SpecCent spCen];
   
   

   %%%% maximum frequency amp
   L = length(zd); NFFT = 2^nextpow2(L);  
   Y = fft(zd,NFFT)/L;
   f = 200/2*linspace(0,1,NFFT/2+1);
   % Plot single-sided amplitude spectrum.      plot(f,2*abs(Y(1:NFFT/2+1)))
   
   [llZ llindZ]=max(2*abs(Y(1:NFFT/2+1)));vvZ = llZ./llindZ; 
   rmsZ = rms(2*abs(Y(1:NFFT/2+1)));
   fW = (abs(Y(1:250)));
   sk1 = (sum(fW-mean(fW))^3)/250;
   sk2 = (sqrt((sum(fW-mean(fW)))^2/250))^3;
   skw = (sk1/sk2);   
   
   L = length(xd); NFFT = 2^nextpow2(L);  
   Y = fft(xd,NFFT)/L; Y(1:3)=0; f = 200/2*linspace(0,1,NFFT/2+1);
   [llX llindX]=max(2*abs(Y(1:NFFT/2+1)));vvX = llX./llindX; 
   rmsX = rms(2*abs(Y(1:NFFT/2+1)));
   
   L = length(yd); NFFT = 2^nextpow2(L);  
   Y = fft(yd,NFFT)/L; Y(1:3)=0; f = 200/2*linspace(0,1,NFFT/2+1);
   [llY llindY]=max(2*abs(Y(1:NFFT/2+1)));vvY = llX./llindX; 
   rmsY = rms(2*abs(Y(1:NFFT/2+1)));
  
   Av = (rmsZ+rmsX+rmsY)/3;
   PFA = (llZ+llX+llY)/3;
   FA = (llindZ+llindX+llindY)/3;
   PF_FA = (vvZ+vvX+vvY)/3;
   
   rmsA = [rmsA Av];            % RMS of frequency amplitude 
   maxPFA = [maxPFA PFA];       % Maximum power of frequency amplitude
   maxFA = [maxFA FA];          % Dominent frequency
   maxPF_FA = [maxPF_FA PF_FA]; % Maximum power of frequency amplitude / dominent frequency
   skwnss = [skwnss skw];    % skewness
   

    %%% Energy
    EZ = feature_energy(zd);
    EX = feature_energy(xd);
    EY = feature_energy(yd);
    E = (EZ + EX + EY)/3;  
    Energy = [Energy E]; % Energy 
 
   
   %%%% polarization
   r=[];
   for i=1:16;
       a=((i-1)*50)+1;
       x=xd(a:i*50);y=yd(a:i*50); z=zd(a:i*50);
       ns=min([length(x) length(y) length(z)]);
       rectilinP=zeros(ns,1);
       MP=cov([x(:) y(:) z(:)]); % covariance(xP,yP,zP);
       [pP,DP] = eig(MP);
       rectilinP =1-((DP(1,1)+DP(2,2))/(2*DP(3,3)));
       r=[r rectilinP];
   end
   [m idx]=max(r);    
   i=idx;a=((i-1)*50)+1;
   x=xd(a:i*50);y=yd(a:i*50); z=zd(a:i*50);
   ns=min([length(x) length(y) length(z)]);
   rectilinP=zeros(ns,1);
   azimuthP=zeros(ns,1);
   dipP=zeros(ns,1);
   MP=cov([x(:) y(:) z(:)]); % covariance(xP,yP,zP);
    [pP,DP] = eig(MP);
        %%% DP contains the eigenvalues of the covariance matrix, with
        %%% DP(3,3)>DP(2,2)>DP(1,1)
        %%% pP contains the eigenvectors, where the first column is the
        %%% eigenvector that corresponds to the smallest eigenvalue, the
        %%% second one to the intermedian eigenvalue and the third one to
        %%% the largest eigenvalue (this one shows the dominant particle motion)
    rectilinP =1-((DP(1,1)+DP(2,2))/(2*DP(3,3)));  
        %%%% the degree of rectilinearity Rec 0 for circular polarization ?1  ?2  ?3,
        %%%% and Rec=1 for rectilinear polarization (?1=1 and ?2=?3=0). P and S waves
        %%%% are body waves, so their degree of rectilinearity is close to 1.
    
%     azimuthP = atan(pP(2,3)/pP(1,3))*180/pi;
      az = atan(pP(2,3)*sign(pP(3,3))/pP(1,3)*sign(pP(3,3)))*180/pi;  % M
      phi = acos(abs(pP(3,3)))*180/pi;   %  M 
      
      
    dipP = atan(pP(3,3)/sqrt(pP(2,3)^2+pP(1,3)^2))*180/pi;% The dip of maximum polarization
    dip_rect=sin(abs(dipP*pi/180)).*rectilinP;
    dipp=sin(abs(dipP*pi/180));
   
    
    rec = [rec rectilinP];
    azimuth = [azimuth az];
    via = [via phi];  % vertical incidence angle 
    DIP = [DIP dipP];
    dd_re = [dd_re dip_rect];
    DIPp = [DIPp dipp];
    
    
    %%% cross correlation
    slave = zd;
%     npvdD=length(masterD);
%     npvdS=length(masterS);
%     npvs=length(slave); 
% 
%     delete dat_s1
%     delete dat_s2
    
%     np=min(npvdS,npvdD);np=min(np,npvs); %Set length of arrays to the smallest 
%     dat_s1(1:np)=masterD(1:np);
%     dat_s2(1:np)=masterS(1:np);
    cross_cor=xcorr(masterD,slave,'coeff');
    [ccor_max, indx]=max(cross_cor);
    ccD=[ccD ccor_max];

    cross_cor=xcorr(masterS,slave,'coeff');
    [ccor_max, indx]=max(cross_cor);
    ccS=[ccS ccor_max];
    
    

%     %%% wavelet Transforming
%     datat=linspace(0,(dt*length(yd)),length(yd));
%     [wl,wlas,wldWx] = cwt_fw(yd,opt.type,opt.nv,data.dt);
%     
%     [na n] = size(wl); col = 0;
%     for i = 1:n
%     col = col + abs(wl(:,i));
%     end
% 
%     xotsu = otsu(na,col); % The Otsu method for finding the optimum point of separation 
%     otY = [otY xotsu];
%     
%     datat=linspace(0,(dt*length(xd)),length(xd));
%     [wl,wlas,wldWx] = cwt_fw(xd,opt.type,opt.nv,data.dt);
%     
%     [na n] = size(wl); col = 0;
%     for i = 1:n
%     col = col + abs(wl(:,i));
%     end
% 
%     xotsu = otsu(na,col); % The Otsu method for finding the optimum point of separation 
%     otX = [otX xotsu];
%     
%     datat=linspace(0,(dt*length(zd)),length(zd));
%     [wl,wlas,wldWx] = cwt_fw(zd,opt.type,opt.nv,data.dt);
%     
%     [na n] = size(wl); col = 0;
%     for i = 1:n
%     col = col + abs(wl(:,i));
%     end
% 
%     xotsu = otsu(na,col); % The Otsu method for finding the optimum point of separation 
%     otZ = [otZ xotsu];
%     
%     data.x=masterD;
%     [wlD,wlasD,wldWxD] = cwt_fw(data.x,opt.type,opt.nv,data.dt);
% 
%     data.x=masterS;
%     [wlS,wlasS,wldWxS] = cwt_fw(data.x,opt.type,opt.nv,data.dt);
%    
%     cD = normxcorr2(abs(wl),abs(wlD));
%     cS = normxcorr2(abs(wl),abs(wlS));
%     
%     xD= mean(mean(abs(cD)));
%     xS= mean(mean(abs(cS)));
%     
%     ccnAb2D=[ccnAb2D xD];
%     ccnAb2S=[ccnAb2S xS];
% 
%     [ypeakD, xpeakD] = find(cD==max(cD(:)));
%     [ypeakS, xpeakS] = find(cS==max(cS(:)));
% 
%     ypicD=[ypicD ypeakD];
%     ypicS=[ypicS ypeakS];
%     
%     cD = normxcorr2(real(wl),real(wlD));
%     cS = normxcorr2(real(wl),real(wlS));
%     
%     xD= mean(mean(abs(cD)));
%     xS= mean(mean(abs(cS)));
%     
%     ccnRel2D=[ccnRel2D xD];
%     ccnRel2S=[ccnRel2S xS];
%     
%     RD = corr2(abs(wl),abs(wlD));
%     RS = corr2(abs(wl),abs(wlS));
%     ccAb2D=[ccAb2D RD];
%     ccAb2S=[ccAb2S RS];
%     
%     RD = corr2(real(wl),real(wlD));
%     RS = corr2(real(wl),real(wlS));
%     ccRel2D=[ccRel2D RD];
%     ccRel2S=[ccRel2S RS];

     
%     %%% discrit wavelet 
%       lev   = 5; wname = 'sym2'; nbcol = 64;
%       [c,l] = wavedec(zd,lev,wname);
%       len = length(zd);
%       cf = zeros(lev,len);
%       for k = 1:lev
%          d = detcoef(c,l,k);
%          d = d(:)';
%          d = d(ones(1,2^k),:);
%          cf(k,:) = wkeep1(d(:)',len);
%       end
%       cf =  cf(:);
%       I = find(abs(cf)<sqrt(eps));
%       cf(I) = zeros(size(I));
%       cf = reshape(cf,lev,len);
% 
%      [c,l] = wavedec(masterD,lev,wname);
%       signal=masterD;
%       len = length(signal);
%       cfdD = zeros(lev,len);
%       for k = 1:lev
%          d = detcoef(c,l,k);
%          d = d(:)';
%          d = d(ones(1,2^k),:);
%          cfdD(k,:) = wkeep1(d(:)',len);
%       end
%       cfdD =  cfdD(:);
%       I = find(abs(cfdD)<sqrt(eps));
%       cfdD(I) = zeros(size(I));
%       cfdD = reshape(cfdD,lev,len);
%       
%       [c,l] = wavedec(masterS,lev,wname);
%       signal=masterS;
%       len = length(signal);
%       cfdS = zeros(lev,len);
%       for k = 1:lev
%          d = detcoef(c,l,k);
%          d = d(:)';
%          d = d(ones(1,2^k),:);
%          cfdS(k,:) = wkeep1(d(:)',len);
%       end
%       cfdS =  cfdS(:);
%       I = find(abs(cfdS)<sqrt(eps));
%       cfd(I) = zeros(size(I));
%       cfdS = reshape(cfdS,lev,len);
%       
%     cxD = normxcorr2(cf,cfdD);
%     cxS = normxcorr2(cf,cfdS);
%     
%     xxD= mean(mean(abs(cxD)));
%     xxS= mean(mean(abs(cxS)));
%     
%     xppD=[xppD xxD];
%     xppS=[xppS xxS];
%     
%     cxD = corr2(real(cf),real(cfdD));
%     cxS = corr2(real(cf),real(cfdS));
% 
%     xpp2D=[xpp2D cxD];
%     xpp2S=[xpp2S cxS];
    

    
%     %%% Coherency 
%     [sD dD]=semblance(datat,slave,masterD,150);
%     [sS dS]=semblance(datat,slave,masterS,150);
%     
%     q = sD > 0.8;qq=sum(sum(q));
%     [a b] = size(sD); n=a*b; sD=qq./n;
%     semSD=[semSD sD];
%     
%     q = sS >= 0.8;qq=sum(sum(q));
%     [a b] = size(sS); n=a*b; sS=qq./n;
%     semSS=[semSS sS];
%     
%     q = dD > 0.75;qq=sum(sum(q));
%     [a b] = size(dD); n=a*b; dD=qq./n;
%     semDD=[semDD sD];
%     
%     q = dS > 0.75;qq=sum(sum(q));
%     [a b] = size(dS); n=a*b; dS=qq./n;
%     semDS=[semDS dS];
%     
%     
%     
%     %%%% Semblence 
%      sSlave = specgram(slave,200);
%      sMasterD = specgram(masterD,200);
%      sMasterS = specgram(masterS,200);
%      semb = [];
%      for i = 2: 100;
%          s1 = ((real(sSlave(i))*real(sMasterD(i)))+ (imag(sSlave(i))*imag(sMasterD(i))));
%          s2 = sqrt((real(sSlave(i))^2*imag(sSlave(i))^2))+ sqrt((real(sMasterD(i))^2*imag(sMasterD(i))^2));
%          sem = s1./s2;
%          semb = [semb sem];
%                   
%      end
%      q = abs(semb) > 0.75;semblD=sum(sum(q));
%      
%      semblanceD = [semblanceD semblD];
%      
%           semb = [];
%      for i = 2: 100;
%          s1 = ((real(sSlave(i))*real(sMasterS(i)))+ (imag(sSlave(i))*imag(sMasterS(i))));
%          s2 = sqrt((real(sSlave(i))^2*imag(sSlave(i))^2))+ sqrt((real(sMasterS(i))^2*imag(sMasterS(i))^2));
%          sem = s1./s2;
%          semb = [semb sem];
%                   
%      end
%      q = abs(semb) > 0.75;semblS=sum(sum(q));
%      semblanceS = [semblanceS semblS];
%      
%      
%      %%%% Envelope 
%      aa = (slave(:)).^2;
%      bb = hilbert(slave(:)).^2;
%      eSlave = sqrt(aa + bb);
%      
%      aa = (masterS(:)).^2;
%      bb = hilbert(masterS(:)).^2;
%      eMasterS = sqrt(aa + bb);
%      
%      aa = (masterD(:)).^2;
%      bb = hilbert(masterD(:)).^2;
%      eMasterD = sqrt(aa + bb);
%      
%      envD =  abs(sum(abs(eSlave-eMasterD))/sum(eSlave));
%      envelopD = [envelopD envD];
%      
%      envS =  abs(sum(abs(eSlave-eMasterS))/sum(eSlave));
%      envelopS = [envelopS envS];
%     
%      
%      
%      %% Euclidian Distance 
%      
%      [fy,f,PDS]=FFT(slave,200);
%      pdSlv = PDS(1:500);
%      
%      [fy,f,PDS]=FFT(masterD,200);
%      pdMstD = PDS(1:500);
%      
%      [fy,f,PDS]=FFT(masterS,200);
%      pdMstS = PDS(1:500);
%      
% 
%      dD = sum((eSlave-eMasterD).^2).^0.5; duD = [duD abs(dD)];
%      dS = sum((eSlave-eMasterS).^2).^0.5; duS = [duS abs(dS)];
%      
%      d2D = (eSlave-eMasterD)'*(eSlave-eMasterD); d2uD = [d2uD d2D];
%      d2S = (eSlave-eMasterS)'*(eSlave-eMasterS); d2uS = [d2uS d2S];
%      

 end 
 % 1-maxAmpZ  2-maxAmpX  3-maxAmpY 4-slopeZ 5-slopeX 6-slopeY 
 % 7-spectralCentroidZ 8-spectralCentroidX 9-spectralCentroidY 
 % 10-maxFreqAmpZ 11-maxFreqAmpX 12-maxFreqAmpY
 % 13-rectilinearity 14-dipp 15-dipp*rec 
 % 16-crossCorrDeep 17-crossCorrShal 18-otsuWtZ 19-otsuWtX 20-otsuWtY
 % 21-2dNormCrossCorAbDeep 22-2dNormCrossCorAbShal 
 % 23-2dNormCrossCorAbDeepYpick 24-2dNormCrossCorAbShalYpic
 % 25-2dNormCrossCorRelDeep 26-2dNormCrossCorRelShal 
 % 27-2dCrossCorAbDeep 28-2dCrossCorAbShal 
 % 29-2dCrossCorRelDeep 30-2dCrossCorRelShal 
 % 31-2dNormCrossCorAbDeepDWT 32-2dNormCrossCorAbShalDWT 
 % 33-2dCrossCorAbDeepRelDWT 34-2dCrossCorRelShalDWT 
 % 35-semblenceSdeep 36-semblenceSshal 37-semblenceDdeep 38-semblenceDshal  
 
%  ev = [maZ' maX' maY' ssz' ssx' ssy' scz' scx' scy'... 
%       mfaZ' mfaX' mfaY' rec' dip' dd_re' ccD' ccS' otZ' otX' otY' ... 
%       ccnAb2D' ccnAb2S' ypicD' ypicS' ccnRel2D' ccnRel2S' ...
%       ccAb2D' ccAb2S' ccRel2D' ccRel2S' xppD' xppS' xpp2D' xpp2S' ...
%       semSD' semSS' semDD' semDS'];
%  
% 
%  
% save('OUT2.mat','ev'); csvwrite('OUT2.csv',ev);




    

% % GISMO PLOT
% 
% ww=[]
% for ii=1:1033
%     disp(ii)
% 
% k = sprintf('%s_%s','BY', num2str(ii));
% 
% v = sprintf('/%s/%s/%s','deepvv2', num2str(ii),'sac');
% [t, d, hdr] = fget_sac(sprintf(v));
%    
% load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'zd.mat'));
% 
% W = waveform;
% sta = num2str(ii);
% cha = num2str(hdr.station.kstnm);
% net = '0';
% loc = '0';
% scnl = scnlobject(strtrim(sta),strtrim(cha),strtrim(net),strtrim(loc));
% W = set(W,'scnlobject',scnl);
% 
% 
% W = set(W,'data',zd);
% W = set(W,'start',0);
% W = set(W,'freq',200);
% 
% ww=[ww W];
% end
% 
% c = correlation(ww)
% c = taper(c);
% c = xcorr(c,[0 3.5]);
%  Hierarchical cluster tree 
%  c = linkage(c);
% plot(c,'den');
% close(gcf)
% set(gcf,'Position',[50 50 600 400]);
