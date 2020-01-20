
clear all
close all

%%%% Defining master events    %%%%%%%%%%%%%%%%%%%%%%%%%%%
 load('masterD.mat'); masterD=d;clear d;dt=1/200; % 74 507
 load('masterS.mat'); masterS=z;clear z;
 
%%% Parameters for Calculate the wavelet transform -
 opt.type = 'morlet';         % Mother wavelet type
 opt.padtype = 'symmetric';   % padded via symmetrization
 opt.rpadded = 1;
 opt.nv = 16; 
 data.dt=1/200;
 dt=1/200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 fileID = fopen('outFeature.csv','w'); 

fprintf(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'...
    ,'evntName','lat','lon','depth','mag','maxAmp','spCen','rmsA','maxPFA' ...
    ,'maxFA','maxPF_FA','Energy','rect','azmth','indAngl','Dip','DipRec','ccD','ccS','xotsu','ccnAb2D' ...
    ,'ccnAb2S','ypicD','ypicS','ccnRel2D','ccnRel2S','ccAb2D','ccAb2S','ccRel2D','ccRel2S' ...
    ,'xp2D','xp2S','xpp2D','xpp2S','semD','semS','skwnss','semblanceD'...
    ,'semblanceS','envelopD','envelopS','udD','udS','ud2D','ud2S','class');

ddd=[]; pfadd=[];spdd=[]; dccdd=[];rectdd=[];dipdd=[];dp=[];mg=[];am=[];spC=[];
    
%  ev = [events(ii).name  lat lon dpth maxAmp spCen rmsA maxPFA maxFA maxPF_FA Energy rect ...
%     dip dip_rect azth ccD ccS xotsu ccnAb2D ccnAb2S ypicD ypicS ccnRel2D ccnRel2S...
%     ccAb2D ccAb2S ccRel2D ccRel2S xp2D xp2S xpp2D xpp2S semD semS class] 
% 
% tline = textread('outFeature.csv','%s','delimiter','\n','whitespace','');
% csvwrite('outFeature.csv',tline);

tline = textread('evntlist.txt','%s','delimiter','\n','whitespace','');

all_events = dir('DATA/*');
ev = all_events([all_events(:).isdir]);
events = all_events(4:length(all_events));
num_dir = numel(events);

 compnm = 0;
 
  for ii=1:num_dir 
     
    disp(events(ii).name)
    clear hdr
    clear d
    v = sprintf('%s/%s/%s','DATA', num2str(events(ii).name),'LA*');
    all_seismograms = dir(sprintf(v));
       
    
    %%% header's info %%%%%%%%%%%%%%%
    vv = sprintf('%s/%s/%s.*','DATA', num2str(events(ii).name),'addev');
    macnm = dir(sprintf(vv));
    tl = textread(macnm.name,'%s','delimiter','\n','whitespace',' ');
    headInfo =strread(char(tl(2)),'%s','delimiter',' ');
    dpth = str2num(char(headInfo(7)));
    lat = str2num(char(headInfo(3)));
    lon = str2num(char(headInfo(5)));
    mag = str2num(char(headInfo(9)));
    ddd=[ddd dpth];
    
    if dpth >= 1000
    class = 'Deep';
    else
    class = 'Shallow';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    stlst=[];
    for jj=1:length(all_seismograms )
    disp(all_seismograms (jj).name)
    seis = all_seismograms (jj).name;
    C = strsplit(seis,'.');
    
    if length(seis) == 11
    stn = sprintf('%s.%s',cell2mat(C(1)),'00');
    stlst=[stlst; stn];
    else 
    stn = sprintf('%s.%s',cell2mat(C(1)),cell2mat(C(4)));
    stlst=[stlst; stn];
    end
    end
     
    stnu = unique(stlst,'rows');
    stnus = cellstr(stnu);
     
    for kk = 1:length(stnus)
    y=strread(char(stnus(kk)),'%s','delimiter','.');
     
    if  str2num(char(y(2))) == 00
    v = sprintf('%s/%s/%s.*Z.*','DATA', num2str(events(ii).name),char(y(1)));
    comn = dir(sprintf(v));
    [t, d, head] = fget_sac(sprintf(comn.name));
    eval(sprintf('%s.Z = d', char(y(1)))); 
    
    v = sprintf('%s/%s/%s.*1.*','DATA', num2str(events(ii).name),char(y(1)));
    comn = dir(sprintf(v));
    [t, d, hdr] = fget_sac(sprintf(comn.name));
    eval(sprintf('%s.X = d', char(y(1)))); 
    
    v = sprintf('%s/%s/%s.*2.*','DATA', num2str(events(ii).name),char(y(1)));
    comn = dir(sprintf(v));
    [t, d, hdr] = fget_sac(sprintf(comn.name));
    eval(sprintf('%s.Y = d', char(y(1))));
    
    else
    v = sprintf('%s/%s/%s.*Z.*.%s','DATA', num2str(events(ii).name),char(y(1)),char(y(2)));
    comn = dir(sprintf(v));
    [t, d, hdr] = fget_sac(sprintf(comn.name));
    eval(sprintf('%s.Z = d', char(y(1)))); 
    
    v = sprintf('%s/%s/%s.*1.*.%s','DATA', num2str(events(ii).name),char(y(1)),char(y(2)));
    comn = dir(sprintf(v));
    [t, d, hdr] = fget_sac(sprintf(comn.name));
    eval(sprintf('%s.X = d', char(y(1)))); 
    
    v = sprintf('%s/%s/%s.*2.*.%s','DATA', num2str(events(ii).name),char(y(1)),char(y(2)));
    comn = dir(sprintf(v));
    [t, d, hdr] = fget_sac(sprintf(comn.name));
    eval(sprintf('%s.Y = d', char(y(1))));
    
    end
    end
    
  
 

    %%% maximum amplitude  %%%%%%%%%%%%%%%%%%%%%%%%%%
    p=0;
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    xd = eval(sprintf('%s.X', char(y(1))));
    yd = eval(sprintf('%s.Y', char(y(1))));   
    
    llz = abs(zd);[llz idz]= max(llz);
    llx = abs(xd);[llx idx]= max(llx);
    lly = abs(yd);[lly idy]= max(lly);
    ll = (llz + llx + lly)/3;
  
   
    p = p + ll;
    end
    maxAmp = p/length(stnus); % max amp

    
    
    
   %%% Spectral Centroid  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    p=0;
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    xd = eval(sprintf('%s.X', char(y(1))));
    yd = eval(sprintf('%s.Y', char(y(1))));   
    
    zd = zd(1:2400);
    xd = xd(1:2400);
    yd = yd(1:2400);
    
   [CZ, CMean, CSD, CMax] = SpecCentroid(zd,200);
   [CX, CMean, CSD, CMax] = SpecCentroid(xd,200);
   [CY, CMean, CSD, CMax] = SpecCentroid(yd,200);
   sp = (CZ + CX + CY)/3;
   
   compnm = compnm +1;
   p = p + sp;
   end
   spCen = p/length(stnus); 
   spdd=[spdd spCen];

    

   %%% maximum frequency amp %%%%%%%%%%%%%%%%%%%%%%%%
    p1=0; p2=0; p3=0; p4=0; p5=0;
    
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    xd = eval(sprintf('%s.X', char(y(1))));
    yd = eval(sprintf('%s.Y', char(y(1))));   
    
   L = length(zd); NFFT = 2^nextpow2(L);  
   Y = fft(zd,NFFT)/L;
   f = 200/2*linspace(0,1,NFFT/2+1);
   [llZ llindZ]=max(2*abs(Y(1:NFFT/2+1)));vvZ = llZ./llindZ; 
   rmsZ = rms(2*abs(Y(1:NFFT/2+1)));
   
   fW = (abs(Y(1:250)));
   sk1 = (sum(fW-mean(fW))^3)/250;
   sk2 = (sqrt((sum(fW-mean(fW)))^2/250))^3;
   skwZ = (sk1/sk2);  
   
   L = length(xd); NFFT = 2^nextpow2(L);  
   Y = fft(xd,NFFT)/L; Y(1:3)=0; f = 200/2*linspace(0,1,NFFT/2+1);
   [llX llindX]=max(2*abs(Y(1:NFFT/2+1)));vvX = llX./llindX; 
   rmsX = rms(2*abs(Y(1:NFFT/2+1)));
   
   fW = (abs(Y(1:250)));
   sk1 = (sum(fW-mean(fW))^3)/250;
   sk2 = (sqrt((sum(fW-mean(fW)))^2/250))^3;
   skwX = (sk1/sk2); 
   
   L = length(yd); NFFT = 2^nextpow2(L);  
   Y = fft(yd,NFFT)/L; Y(1:3)=0; f = 200/2*linspace(0,1,NFFT/2+1);
   [llY llindY]=max(2*abs(Y(1:NFFT/2+1)));vvY = llX./llindX; 
   rmsY = rms(2*abs(Y(1:NFFT/2+1)));
   
   fW = (abs(Y(1:250)));
   sk1 = (sum(fW-mean(fW))^3)/250;
   sk2 = (sqrt((sum(fW-mean(fW)))^2/250))^3;
   skwY = (sk1/sk2); 
  
   Av = (rmsZ+rmsX+rmsY)/3;
   PFA = (llZ+llX+llY)/3;
   FA = (llindZ+llindX+llindY)/3;
   PF_FA = (vvZ+vvX+vvY)/3;
   skw = (skwZ +skwX +skwY )/3;
   
   p1 = p1 +  Av;
   p2 = p2 +  PFA;
   p3 = p3 +  FA;
   p4 = p4 +  PF_FA;
   p5 = p5 + skw;
    end
    
   rmsA = p1/length(stnus);        % RMS of frequency amplitude
   maxPFA = p2/length(stnus);      % Maximum power of frequency amplitude
   maxFA = p3/length(stnus);       % Dominent frequency
   maxPF_FA = p4/length(stnus);    % Maximum power of frequency amplitude / dominent frequency
   skwnss = p5/length(stnus);      % average power spectral skewness around dominant frequency

   pfadd = [pfadd maxPFA];
   
   %%% Energy Density   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p=0; 
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    xd = eval(sprintf('%s.X', char(y(1))));
    yd = eval(sprintf('%s.Y', char(y(1))));  
    
    zd = zd(1:2400);
    xd = xd(1:2400);
    yd = yd(1:2400);
    
    EZ = feature_energy(zd);
    EX = feature_energy(xd);
    EY = feature_energy(yd);
    E = (EZ + EX + EY)/3; 
    p = p + E;
   end
   Energy = p/length(stnus); 
   



   %%% polarization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1=0; p2=0; p3=0; p4=0; p5=0;
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    xd = eval(sprintf('%s.X', char(y(1))));
    yd = eval(sprintf('%s.Y', char(y(1))));  
   
    
    
    % finding the approximate window of P
   r=[];
   for i=1:48;
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
 
   rectilinP =1-((DP(1,1)+DP(2,2))/(2*DP(3,3)));  
   az = atan(pP(2,3)*sign(pP(3,3))/pP(1,3)*sign(pP(3,3)))*180/pi;
   phi = acos(abs(pP(3,3)))*180/pi;
   dipP = atan(pP(3,3)/sqrt(pP(2,3)^2+pP(1,3)^2))*180/pi;% The dip of maximum polarization
   dip_rect = sin(abs(dipP*pi/180)).*rectilinP;
       
   p1 = p1 +  rectilinP;   
   p2 = p2 +  az;       
   p3 = p3 +  phi;
   p4 = p4 +  dipP;
   p5 = p5 +  dip_rect;

 end   
   rect = p1/length(stnus);        % the degree of rectilinearity 
   azmth = p2/length(stnus);         % azimuth
   indAngl = p3/length(stnus);     % incident angle
   Dip = p4/length(stnus);        % Dip
   DipRec = p3/length(stnus);     % Dip * rect

   rectdd= [rectdd rect];
   dipdd = [dipdd Dip]; 

   
   
   %%% cross correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1=0;p2=0;
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    zd = zd(1:2400);
    
    slave = zd;
    npvdD=length(masterD);
    npvdS=length(masterS);
    npvs=length(slave); 

    np=min(npvdS,npvdD);np=min(np,npvs); %Set length of arrays to the smallest 
    dat_s1(1:np)=masterD(1:np);
    dat_s2(1:np)=slave(1:np);
    cross_cor=xcorr(dat_s1,dat_s2,'coeff');
    [ccD, indx]=max(cross_cor);
    p1 = p1 + ccD;

    clear dat_s1;
    dat_s1(1:np)=masterS(1:np);
    cross_cor=xcorr(dat_s1,dat_s2,'coeff');
    [ccS, indx]=max(cross_cor);
    p2 = p2 + ccS;
    
   end
   ccD = p1/length(stnus); 
   ccS = p2/length(stnus); 
    
   




   %%%  CWT wavelet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1 = 0; p3 = 0; p4 = 0;
    p5 = 0; p6 = 0; p7 = 0;
    p8 = 0; p9 = 0; p10 = 0;
    p11 = 0;p12 = 0;
    
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    xd = eval(sprintf('%s.X', char(y(1))));
    yd = eval(sprintf('%s.X', char(y(1))));
    
    zd = zd(1:2400);
    xd = xd(1:2400);
    yd = yd(1:2400);
    
    datat=linspace(0,(dt*length(zd)),length(zd));
    [wl,wlas,wldWx] = cwt_fw(zd,opt.type,opt.nv,data.dt);
    
    [na n] = size(wl); col = 0;
    for i = 1:n
    col = col + abs(wl(:,i));
    end

    xotsu = otsu(na,col); 
    p1 = p1 + xotsu; % xotsu
    
    datat=linspace(0,(dt*length(xd)),length(xd));
    [wl,wlas,wldWx] = cwt_fw(xd,opt.type,opt.nv,data.dt);
    
    [na n] = size(wl); col = 0;
    for i = 1:n
    col = col + abs(wl(:,i));
    end

    xotsu = otsu(na,col); 
    p1 = p1 + xotsu; % xotsu
    
    datat=linspace(0,(dt*length(yd)),length(yd));
    [wl,wlas,wldWx] = cwt_fw(yd,opt.type,opt.nv,data.dt);
    
    [na n] = size(wl); col = 0;
    for i = 1:n
    col = col + abs(wl(:,i));
    end

    xotsu = otsu(na,col); 
    p1 = p1 + xotsu; % xotsu
    
    data.x=masterD;
    [wlD,wlasD,wldWxD] = cwt_fw(data.x,opt.type,opt.nv,data.dt);

    data.x=masterS;
    [wlS,wlasS,wldWxS] = cwt_fw(data.x,opt.type,opt.nv,data.dt);
   
    cD = normxcorr2(abs(wl),abs(wlD));
    cS = normxcorr2(abs(wl),abs(wlS));
    
    xD= mean(mean(abs(cD)));
    xS= mean(mean(abs(cS)));
    p3 = p3 + xD;  %  ccnAb2D
    p4 = p4 + xS;  %  ccnAb2S

    [ypeakD, xpeakD] = find(cD==max(cD(:)));
    [ypeakS, xpeakS] = find(cS==max(cS(:)));
    p5 = p5 + ypeakD; % ypicD
    p6 = p6 + ypeakS; % ypicS
    
    cD = normxcorr2(real(wl),real(wlD));
    cS = normxcorr2(real(wl),real(wlS));
    
    xD = mean(mean(abs(cD)));
    xS = mean(mean(abs(cS)));
    p7 = p7 + xD; % ccnRel2D
    p8 = p8 + xS; % ccnRel2S
    
    RD = corr2(abs(wl),abs(wlD));
    RS = corr2(abs(wl),abs(wlS));
    p9 = p9 + RD; % ccAb2D
    p10 = p10 + RS; % ccAb2S
    
    RD = corr2(real(wl),real(wlD));
    RS = corr2(real(wl),real(wlS));
    p11 = p11 + RD; % ccRel2D
    p12 = p12 + RS; % ccRel2S
    
    end
    
    xotsu = p1/3*length(stnus);   % xotsu
    ccnAb2D = p3/length(stnus);   %  ccnAb2D
    ccnAb2S = p4/length(stnus);   %  ccnAb2S
    ypicD = p5/length(stnus);   % ypicD
    ypicS = p6/length(stnus);   % ypicS
    ccnRel2D = p7/length(stnus);   % ccnRel2D
    ccnRel2S = p8/length(stnus);   % ccnRel2S
    ccAb2D = p9/length(stnus);   % ccAb2D
    ccAb2S = p10/length(stnus); % ccAb2S
    ccRel2D = p11/length(stnus); % ccRel2D
    ccRel2S = p12/length(stnus); % ccRel2S
    
     dccdd = [dccdd ccnAb2D];
    
    %%% discrit wavele %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1 = 0; p2 = 0; p3 = 0; p4=0;
    
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    zd = zd(1:2400);

    
    lev   = 5; wname = 'sym2'; nbcol = 64;
    [c,l] = wavedec(zd,lev,wname);
    len = length(zd);
    cf = zeros(lev,len);
    for k = 1:lev
        d = detcoef(c,l,k);
        d = d(:)';
        d = d(ones(1,2^k),:);
        cf(k,:) = wkeep1(d(:)',len);
    end
    
    cf =  cf(:);
    I = find(abs(cf)<sqrt(eps));
    cf(I) = zeros(size(I));
    cf = reshape(cf,lev,len);

    [c,l] = wavedec(masterD,lev,wname);
    signal=masterD;
    len = length(signal);
    cfdD = zeros(lev,len);
    for k = 1:lev
       d = detcoef(c,l,k);
       d = d(:)';
       d = d(ones(1,2^k),:);
       cfdD(k,:) = wkeep1(d(:)',len);
    end
      cfdD =  cfdD(:);
      I = find(abs(cfdD)<sqrt(eps));
      cfdD(I) = zeros(size(I));
      cfdD = reshape(cfdD,lev,len);
      
      [c,l] = wavedec(masterS,lev,wname);
      signal=masterS;
      len = length(signal);
      cfdS = zeros(lev,len);
      for k = 1:lev
         d = detcoef(c,l,k);
         d = d(:)';
         d = d(ones(1,2^k),:);
         cfdS(k,:) = wkeep1(d(:)',len);
      end
      cfdS =  cfdS(:);
      I = find(abs(cfdS)<sqrt(eps));
      cfd(I) = zeros(size(I));
      cfdS = reshape(cfdS,lev,len);
      
    cxD = normxcorr2(cf,cfdD);
    cxS = normxcorr2(cf,cfdS);
    
    xxD= mean(mean(abs(cxD)));
    xxS= mean(mean(abs(cxS)));
    
    p1 = p1 + xxD;
    p2 = p2 + xxS;
    
    cxD = corr2(real(cf),real(cfdD));
    cxS = corr2(real(cf),real(cfdS));

    p3 = p3 + cxD; % xpp2D
    p4 = p4 + cxS; % xpp2S

end 
    xp2D = p1/length(stnus);
    xp2S = p2/length(stnus);
    xpp2D = p3/length(stnus);
    xpp2S = p4/length(stnus);

    
    %%% spectral coherancy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1 = 0; p2 = 0;
    
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    zd = zd(1:2400);
    
    [sD dD]=semblance(datat,zd(1:np),masterD(1:np),150);
    [sS dS]=semblance(datat,zd(1:np),masterS(1:np),150);
    
    q = sD > 0.8;qq=sum(sum(q));
    [a b] = size(sD); n=a*b; sD=qq./n;
    p1 = p1 + sD;
    
    q = sS >= 0.8;qq=sum(sum(q));
    [a b] = size(sS); n=a*b; sS=qq./n;
    p2 = p2 + sS;

 end 
    semD = p1/length(stnus);
    semS = p2/length(stnus);

    
    %%% spectral semblance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1 = 0; p2 = 0;
    
    sMasterD = specgram(masterD,200);
    sMasterS = specgram(masterS,200);
     
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    zd = zd(1:2400);
    
    sSlave = specgram(zd,200);

    
    semb = [];
     for i = 2: 100;
         s1 = ((real(sSlave(i))*real(sMasterD(i)))+ (imag(sSlave(i))*imag(sMasterD(i))));
         s2 = sqrt((real(sSlave(i))^2*imag(sSlave(i))^2))+ sqrt((real(sMasterD(i))^2*imag(sMasterD(i))^2));
         sem = s1./s2;
         semb = [semb sem];
                  
     end
    q = abs(semb) > 0.75;semblD=sum(sum(q));
    p1 = p1 + semblD;
    
    semb = [];
     for i = 2: 100;
         s1 = ((real(sSlave(i))*real(sMasterS(i)))+ (imag(sSlave(i))*imag(sMasterS(i))));
         s2 = sqrt((real(sSlave(i))^2*imag(sSlave(i))^2))+ sqrt((real(sMasterS(i))^2*imag(sMasterS(i))^2));
         sem = s1./s2;
         semb = [semb sem];
                  
     end
     q = abs(semb) > 0.75;semblS=sum(sum(q));
     p2 = p2 + semblS;

 end 
    semblanceD = p1/length(stnus);
    semblanceS = p2/length(stnus);
    
    
    
    
    %%%% Envelope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1 = 0; p2 = 0;
    
     aa = (masterS(:)).^2;
     bb = hilbert(masterS(:)).^2;
     eMasterS = sqrt(aa + bb);
     
     aa = (masterD(:)).^2;
     bb = hilbert(masterD(:)).^2;
     eMasterD = sqrt(aa + bb);
     
     
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    zd = zd(1:2400);
    
    aa = (zd(:)).^2;
    bb = hilbert(zd(:)).^2;
    eSlave = sqrt(aa + bb);

    envD =  abs(sum(abs(eSlave-eMasterD))/sum(eSlave));
    p1 = p1 + envD;
    
    envS =  abs(sum(abs(eSlave-eMasterS))/sum(eSlave));
    p2 = p2 + envS ;

 end 
    envelopD = p1/length(stnus);
    envelopS = p2/length(stnus);
    
    
    
    
     %%%%% spectral distance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1 = 0; p2 = 0;
    
     [fy,f,PDS]=FFT(masterD,200);
     pdMstD = PDS(1:500);
     
     [fy,f,PDS]=FFT(masterS,200);
     pdMstS = PDS(1:500);
    
    for kk = 1:length(stnus)
    y = strread(char(stnus(kk)),'%s','delimiter','.');
    zd = eval(sprintf('%s.Z', char(y(1))));
    zd = zd(1:2400);
    
    [fy,f,PDS]=FFT(zd,200);
    pdSlv = PDS(1:500);
    
    dD = sum((pdSlv-pdMstD).^2).^0.5; 
    p1 = p1 +  abs(dD);
    
    dS = sum((pdSlv-pdMstS).^2).^0.5;
    p2 = p2 + abs(dS);
    
    d2D = (pdSlv-pdMstD)'*(pdSlv-pdMstD);
    p3 = p3 + abs(d2D);
        
    d2S = (pdSlv-pdMstS)'*(pdSlv-pdMstS);
    p4 = p4 + abs(d2S);
    
 end 
    udD = p1/length(stnus)
    udS = p2/length(stnus)
    ud2D = p3/length(stnus)
    ud2S = p4/length(stnus)
    
    
fprintf(fileID,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s \n'...
    , events(ii).name, lat, lon, dpth, mag, maxAmp, spCen, rmsA, maxPFA, maxFA ...
    , maxPF_FA, Energy, rect, azmth, indAngl, Dip, DipRec, ccD, ccS, xotsu ...
    , ccnAb2D, ccnAb2S, ypicD, ypicS, ccnRel2D, ccnRel2S, ccAb2D, ccAb2S...
    , ccRel2D, ccRel2S, xp2D, xp2S, xpp2D, xpp2S, semD, semS, skwnss...
    , semblanceD, semblanceS, envelopD, envelopS, udD, udS, ud2D, ud2S, class);
end 
 fclose(fileID);
 

% %  ev = [maZ' maX' maY' ssz' ssx' ssy' scz' scx' scy'... 
% %       semSD' semSS' semDD' semDS'];

% % save('OUT2.mat','ev'); csvwrite('OUT2.csv',ev);

 
