clc
clear all
close all

% Parameters for Calculate the wavelet transform -
opt.type = 'morlet';         % Mother wavelet type
opt.padtype = 'symmetric';   % padded via symmetrization
opt.rpadded = 1;
opt.nv = 16;  

%load(data.nm)
% data.x=dd;
% data.t=linspace(0,(dt*length(dd)),length(dd));
% data.dt=1/200;
% dt=1/200;


%% SST
data.nm ='2_20150220_040140-17-2';
[data.t,data.x,data.hdr] = fget_sac(data.nm);
data.t = linspace(0,(data.hdr.times.e - data.hdr.times.b),length(data.x));
data.dt = data.t(2)-data.t(1);
% data.x=hipass(data.x,200,2,2,2,'butter');

% wavelet Transforming
[wl.Coef,wl.as,wl.dWx] = cwt_fw(data.x,opt.type,opt.nv,data.dt);
wl.w = imag(wl.dWx ./ wl.Coef / (2*pi));
[ss.Coef,ss.fs] = synsq_cwt_squeeze(wl.Coef, wl.w, data.t, opt.nv, opt);


data.x=data.x/max(data.x);
figure
hold on

subplot 211
plot(data.t,data.x); axis tight

grid on
grid minor
title(data.nm,'Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
xlim([min(data.t) max(data.t)]);
ax = gca;
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';

clear title xlabel ylabel



subplot 212
imagesc(data.t,fliplr(wl.as),abs(flipud(wl.Coef)));

title({'CWT Scalogram'},'Rotation',0,'FontSize',13); 
xlabel({'Time (s)'},'FontSize',11);
ylabel('Scale (1/Hz)','FontSize',11);

set(gca,'YDir','normal');

ax = gca;
ax.CLim = [-1000 10000];
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';

axis tight;
hold off
clear title xlabel ylabel



subplot 412
tplot(wl.Coef,data.t, wl.as); 


subplot 413
tplot(ss.Coef,data.t, ss.fs); 

subplot 414
spectrogram(data.x,30,'yaxis'); colormap (jet)
ax.CLim = [-110 -20];




 




%% detrending
for i = 1:length(d)
nm = d(i).name;
nm=cellstr(nm);
nm = cell2mat(nm);
load(sprintf(nm));
dd=detrend(dd);
save(sprintf(nm),'dd')
end














plot(f,2*abs(Y(1:NFFT/2+1))) 

ax = gca;
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
grid on

%% modified one 
nn='mat';
l=[];
for i = 1:535
load(sprintf('%d.%s%s',i,nn));
L = length(dd);                     % Length of signal
NFFT = 2^nextpow2(L);
Y = fft(dd,NFFT)/L;
Y(1:3)=0;
f = 200/2*linspace(0,1,NFFT/2+1);
[ccor_max, indx]=max(2*abs(Y(1:NFFT/2+1)));
ll=ccor_max./indx;
l=[l ll];
end
l=l';






%% cross correlating 
clear all
load('316.mat')
master=dd;
npvd=length(master);
dt=1/200;

nn='mat';
mcc=[]; zcc=[];
for i = 1:535
load(sprintf('%d.%s%s',i,nn));
slave=dd;

%  Set length of arrays to the smaller of data and synthetic
npvs=length(slave);
np=min(npvd,npvs);

dat_s1(1:np)=master(1:np);
dat_s2(1:np)=slave(1:np);

% tplot=linspace(0,(np-1).*dt,np);
% t=linspace(-dt.*(np),dt.*np,2*np-1);
% nzero=np;

cross_cor=xcorr(dat_s1,dat_s2,'coeff');
[ccor_max, indx]=max(cross_cor);
c_coef=ccor_max;   % maximum of the normalized cross correlation
% c_zero=cross_cor(nzero); % value of the correlation at zero time shift
% tshift=(indx - nzero).*dt; % time shift of the maximum normalized correlation

mcc=[mcc c_coef];
% zcc=[zcc c_zero];


clear npvs
clear np
clear dat_s1
clear dat_s2

end
plot(mcc)


figure;
    subplot(3,1,1);
    plot(tplot,dat_s1,'-k');
    subplot(3,1,2);
    plot(tplot,dat_s2,'-k');
    subplot(3,1,3);
    t1=t(indx);
    plot(t,cross_cor,'-k',t1,c_coef,'ok');
    
  
    
    
%% cross correlating spectrums 
clear all
load('316.mat')
master=dd;
dt=1/200;

L = length(master);
NFFT = 2^nextpow2(L);
Y = fft(master,NFFT)/L;
f =200/2*linspace(0,1,NFFT/2+1);
Ym=2*abs(Y(1:NFFT/2+1));

nn='mat';
mcc=[]; zcc=[];
for i = 1:535
load(sprintf('%d.%s%s',i,nn));
slave=dd;
L = length(slave);
NFFT = 2^nextpow2(L);
Y = fft(slave,NFFT)/L;
f =200/2*linspace(0,1,NFFT/2+1);
Ys=2*abs(Y(1:NFFT/2+1));


    



%% discrit wavelet 
lev   = 5;
wname = 'sym2';
nbcol = 64;
[c,l] = wavedec(data.x,lev,wname);
signal=data.x;
len = length(signal);
cfd = zeros(lev,len);
for k = 1:lev
    d = detcoef(c,l,k);
    d = d(:)';
    d = d(ones(1,2^k),:);
    cfd(k,:) = wkeep1(d(:)',len);
end
cfd =  cfd(:);
I = find(abs(cfd)<sqrt(eps));
cfd(I) = zeros(size(I));
cfd = reshape(cfd,lev,len);
cfd = wcodemat(cfd,nbcol,'row');
colormap(pink(nbcol));
image(cfd);
tics = 1:lev;
labs = int2str((1:lev)');
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YDir = 'normal';
ax.Box = 'On';
ax.YTick = tics;
ax.YTickLabel = labs;
title('Discrete Wavelet Transform, Absolute Coefficients.');
xlabel('Time (or Space)');
ylabel('Level');
title('766-shallow, Absolute Coefficients.');





%% simple bash
j=0
for i in $(ls); do 
    
    cd $i 
        cp LA*.*HZ.*  sac
    cd ..
done

j=0
for i in $(ls); do 
    j=j+1;
    mv -v $i  /cc/$j
 done


 

%% Picking 
for i = 1:length(d)
    disp(d(i).name);
    ii = d(i).name;
    cd(sprintf('%s', ii ));
    
    dd = dir('LA*');
    for k = 1:length(dd)
    iii = dd(k).name;
    C = strsplit(iii,'.');
    
    dz = dir(sprintf('%s.%s.%s.%s',cell2mat(C(1)),'*HZ',cell2mat(C(3)),cell2mat(C(4))));   
    dx = dir(sprintf('%s.%s.%s.%s',cell2mat(C(1)),'*H1',cell2mat(C(3)),cell2mat(C(4))));
    dy = dir(sprintf('%s.%s.%s.%s',cell2mat(C(1)),'*H2',cell2mat(C(3)),cell2mat(C(4))));    

    [data.t, z, data.hdr] = fget_sac(sprintf('%s',dz.name));
    [data.t, x, data.hdr] = fget_sac(sprintf('%s',dx.name));
    [data.t, y, data.hdr] = fget_sac(sprintf('%s',dy.name));

    z = hipass(z,200,2,2,1,'butter','linear');
    x = hipass(x,200,2,2,1,'butter','linear');
    y = hipass(y,200,2,2,1,'butter','linear');

    z = detrend(z);
    x = detrend(x);
    y = detrend(y);
  
    zz = sprintf('%s%s',cell2mat(C(1)),'Z');
    xx = sprintf('%s%s',cell2mat(C(1)),'X');
    yy = sprintf('%s%s',cell2mat(C(1)),'Y');

    Z.(zz) = z;
    X.(xx) = x;
    Y.(yy) = y;
    end
    
    nameZ = fieldnames(Z);
      
    figure (1)
    for m = 1:length(dd)
        subplot(length(dd),1,m)
        mm = nameZ(m);
        plot(data.t,Z.(sprintf('%s',cell2mat(mm))));
    end  
     
end 

    subplot 311; plot(z)
    subplot 312; plot(x)
    subplot 313; plot(y)
    
    prompt = 'Origin? '; o = input(prompt)
    prompt = 'End? '; e = input(prompt)
    
    zz=z(o-200:e);
    xx=x(o-200:e);
    yy=y(o-200:e);   
        
    z=zeros(2000,1); z(1:length(zz))=zz;
    x=zeros(2000,1); x(1:length(xx))=xx;
    y=zeros(2000,1); y(1:length(yy))=yy;
    
% parts = strread(dz.name,'%s','delimiter','.')
    save('z.mat','z')
    save('x.mat','x')
    save('y.mat','y')
    
    cd .. 
    end



    
%% Viewing the mat picked
for i=1018:1033
    disp(i)
    cd(sprintf('%d', i));
    z = load('z.mat');
    x = load('x.mat');
    y = load('y.mat');
    
%     z.dz = hipass(z.dz,200,15,2,1,'butter','linear');
    
    figure (1)
    subplot 311
    plot(z.z)
    
    subplot 312
    plot(x.x)
    
    subplot 313
    plot(y.y)
    
   rr = input('Next?')
 
   cd .. 
end    
    
    



%% Adding zero to the end of tace 
for i=1:225
    disp(i)
    cd(sprintf('%d', i));
    z = load('z.mat');
    x = load('x.mat');
    y = load('y.mat');
    
    dz=zeros(2000,1);
    mi=2000-length(z.z);
    dz(1:length(z.z))=z.z;
    dz(length(z.z)+1:length(z.z)+mi)=0;
    save('z.mat','dz')
    
    dx=zeros(2000,1);
    dx(1:length(x.x))=x.x;
    dx(length(x.x)+1:length(x.x)+mi)=0;
    save('x.mat','dx')
    
    dy=zeros(2000,1);
    dy(1:length(y.y))=y.y;
    dy(length(y.y)+1:length(y.y)+mi)=0;
    save('y.mat','dy')
    
    rr = input('Next?')
   cd .. 
end






%% Manual quality check 
all_files = dir;
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir);
d1 = all_dir(3:num_dir)
badev = [];
fileID = fopen('badev.txt','w');

for i = 1:length(d1)
    disp(i/length(d1));
    disp(d1(i).name);
    ii = d1(i).name;
    cd(sprintf('%s', ii ));
    clear Z
    dd = dir('LA*');
    for k = 1:length(dd)
    iii = dd(k).name;
    C = strsplit(iii,'.');
    [data.t, z, data.hdr] = fget_sac(sprintf('%s',iii));

    z = hipass(z,200,2,2,1,'butter','linear');
    z = detrend(z);
    
    if length(iii) == 11
    zz = sprintf('%s%s%s',cell2mat(C(1)),cell2mat(C(2)),'00');
    else 
    zz = sprintf('%s%s%s',cell2mat(C(1)),cell2mat(C(2)),cell2mat(C(4)));
    end
    
    Z.(zz) = z;
    
    end
    
    nameZ = [];
    nameZ = fieldnames(Z);
      
    figure (1)
    for m = 1:length(dd)
        subplot(length(dd),1,m)
        mm = nameZ(m);
        
        plot(Z.(sprintf('%s',cell2mat(mm))));
        
        ax = gca;
        h = get(gca,'YLabel');
        pos = get(h,'position');
        pos(1) = pos(1) - 0.6 ; 
        set(h,'position',pos)
        
        ylabel(sprintf('%s',cell2mat(mm)),'FontSize',11)
        set(get(gca,'YLabel'),'Rotation',0)
        ax.YAxisLocation = 'right';
       
    end  
        
    r = input('BadStations? 0 or [21 10]')
    
    if(~exist('r'))
        
    for g = length(r)
    delete(sprintf('%s%s.*','LA',num2str(r(g)))); 
    end 
    end
    
    rr = input('Resize? yes=1 No=0')
    if rr == 1
    n = input('Start of Window?')
       figure (1)
    for m = 1:length(dd)
        subplot(length(dd),1,m)
        mm = nameZ(m);
        w = Z.(sprintf('%s',cell2mat(mm)));
        plot(w(n:n+1000));
        
        ax = gca;
        h = get(gca,'YLabel');
        pos = get(h,'position');
        pos(1) = pos(1) - 0.6 ; 
        set(h,'position',pos)
        
        ylabel(sprintf('%s',cell2mat(mm)),'FontSize',11)
        set(get(gca,'YLabel'),'Rotation',0)
        ax.YAxisLocation = 'right';
       
    end 
    end
    
   rrr = input('Do you want to keep it? no=0 yes=1')
    
    if rrr == 0
     S = char(ii)   
     fprintf(fileID,'%s \n',S);
    end
    
   rrrr = input('Next?')
    
    cd .. 
end
fclose(fileID);
system('./deletbad.sh ')



%% denoising

d = dir('untitled');
d1 = d(4:length(d))
for i = 1:length(d1)
    disp(d1(i).name);
    ii = d1(i).name;
    
    v = sprintf('/%s/%s/%s','untitled', num2str(ii),'sac');
    [t, z, hdr] = fget_sac(sprintf(v));
    
    z = load(sprintf('/%s/%s/%s','untitled', num2str(ii),'z.mat'));
    x = load(sprintf('/%s/%s/%s','untitled', num2str(ii),'x.mat'));
    y = load(sprintf('/%s/%s/%s','untitled', num2str(ii),'y.mat'));
    
    zd=denois(z.z,hdr);
    xd=denois(x.x,hdr);
    yd=denois(y.y,hdr);
    
    fz=['untitled/' num2str(ii) '/' 'zd.mat']; 
    fx=['untitled/' num2str(ii) '/' 'xd.mat']; 
    fy=['untitled/' num2str(ii) '/' 'yd.mat']; 
    
    save(fz,'zd');
    save(fx,'xd');
    save(fy,'yd');
      
    
end






%% viewing denoised one and the original
 for ii=1:218
    disp(ii)
    
    v = sprintf('/%s/%s/%s','deepvv2', num2str(ii),'sac');
    [t, z, hdr] = fget_sac(sprintf(v));
    
    z = load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'z.mat'));
    x = load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'x.mat'));
    y = load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'y.mat'));
    
    
    dz = load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'zd.mat'));
    dx = load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'xd.mat'));
    dy = load(sprintf('/%s/%s/%s','deepvv2', num2str(ii),'yd.mat'));
%     z.dz = hipass(z.dz,200,15,2,1,'butter','linear');
    
    figure (1)
    subplot 321
    plot(z.z)
    
    subplot 322
    plot(dz.zd)
    
    subplot 323
    plot(x.x)
    
    subplot 324
    plot(dx.xd)
    
    subplot 325
    plot(y.y)
    
    subplot 326
    plot(dy.yd)
    
    rr = input('Next?')
 
 end  

  
 
 
 
 
%% making ev_list 
 co='bcrmgky';
pink=[1 102/255 1];
purple=[186/255 0 1];
events=dir('2*');
 
 
 fname='export.out';
 tline = textread(fname,'%s','delimiter','\n','whitespace','');
 ne=length(tline);
p=1;
for i=1:ne
y=strread(char(tline(i)),'%s','delimiter',',');
evlat(p)=str2num(char(y(5)));
evlon(p)=str2num(char(y(6)));
evz(p)=str2num(char(y(10)));
evmag(p)=str2num(char(y(11))); 

ss=char(y(1));
yyy=strread(char(y(2)),'%s','delimiter',':');

yy=str2num(ss(1:4));
mm=str2num(ss(6:7));
dd=str2num(ss(9:10));
hh=str2num(char(yyy(1)));
mn=str2num(char(yyy(2)));
rs=str2num(char(yyy(3)));
et(p)=datenum(yy,mm,dd,hh,mn,rs);

if evz(p) >= 1000
%if evz(p) >= 1000 & evmag(p) < +0.1 & evmag(p) > -0.5
p=p+1;
end
end
pp=length(et);


FN=['SaltDome.kml'];
[lat,lon,z] = read_kml(FN);
line(lon,lat,'Color','r','LineWidth',2)
hold on

FN=['sinkhole.kml'];
[lat,lon,z] = read_kml(FN);
line(lon,lat,'Color','b','LineWidth',2)


FN=['sta.dat'];
[sla,slo,sta_n]=textread(FN,'%f %f %s');
scatter(slo,sla,'SizeData',75,'Marker','s','MarkerFaceColor',purple,'MarkerEdgecolor','k');

%PLOT EVENT LOCATION
for i=1:pp;
scatter3(evlon(i),evlat(i),-evz(i),'SizeData',25,'MarkerFaceColor',co(2),'MarkerEdgecolor','k');
end


fid=fopen('ev_list','w');
nzv=length(et)
for i=1:nzv
fprintf(fid,'%s \n',datestr(et,'yyyymmdd_HHMMSS'));
end
fclose(fid);

fid=fopen('area.out','w');

for i=1:nzv
fprintf(fid,'%s %f %f %f \n',datestr(et,'yyyymmdd_HHMMSS'),x(i),y(i),z(i));
end
fclose(fid);
 
 



