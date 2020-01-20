
clc
clear all
close all 
tic

data.nm ='80/LA10.EHZ.YC.02';   % waveform name
% Parameters for Calculate the wavelet transform -
opt.type = 'hshannon';         % Mother wavelet type; 'gauss' 'cmhat'   'morlet'   'hshannon'    'hhhat'    'bump' 
opt.padtype = 'symmetric';   % padded via symmetrization
opt.rpadded = 1;
opt.nv = 64;                 % Number of voices

% Read the original 
h = waitbar(0,'Loading...');
[data.t,data.x,data.hdr] = fget_sac(data.nm);
data.t = linspace(0,(data.hdr.times.e - data.hdr.times.b),length(data.x));
data.dt = data.t(2)-data.t(1);
close(h)
clear h

%%  wavelet Transforming

h = waitbar(0.1,'Wavelet Transforming...');

[wl.Coef,wl.as,wl.dWx] = cwt_fw(data.x,opt.type,opt.nv,data.dt);

waitbar(0.8,h,'Wavelet Transforming...'); close(h)

%% Segmentation 
[dec nn ]= segment(wl.Coef, wl.dWx,wl.as,data);  

%% Synchronize Denoising 
[denoised ] = synchDenois(dec,opt,data,nn);

%% Post Denoising  
% forward wavelet transform
[wdn,w.as] = cwt_fw(denoised,opt.type,opt.nv,data.dt);

% Customized thresholding 
[ dnCuCoef ] = customThrSoft(wdn,opt,nn);
dnCuCoef(isnan(dnCuCoef)) = 0;
x = cwt_iw(dnCuCoef, opt.type, opt, opt.nv); 

% %% Detection
% [TT, ff, aa] = synsq_cwt_fw(data.t, x, opt.nv , opt.type);
% 
% % characteristic function DF
% [na n] = size(TT); 
% ee  = zeros(na,n);
% for i = 1:na
% v = real(TT(i,:));
% a = (v).^2 ;
% b = (hilbert(v)).^2;
% ee(i,:) = sqrt(a+b);
% end
% row = sum(abs(ee));
% 
% rr = row';
% L = 30;
% z =zeros(L,1);
% rowm =[z; rr; z];
% er = [];
% 
%  for i = 1: length(rowm)-2*L;
%      p = L+i;
%      rowUp = rowm(p:p+L);
%      rowUp = rowUp.^2;
%      
%      rowDn = rowm(p-L:p);
%      rowDn = rowDn.^2;
%      
%      e= sum(rowUp)/sum(rowDn);
%      
%      er = [er;e];
% 
%  end
%  
% [pks,locs] = findpeaks(er);
% out = []; T = 0.15*max(er);
% for f=1:length(pks);
%     if pks(f)> T
%        out = [out locs(f)];
%     end
% end
% 
% up.T = zeros(size(out)); 
% for i = 1:length(out);
%     cut = er(out(i)-30:out(i));
%     DV=movingslope(cut,2,1,data.dt);
%     [m idx] = min(abs(DV));
%     up.T(i) = out(i)-(30-i);
%  
% end



% xr = smooth(row,0.01,'loess');
% g=f_cumul(xr);
% 
% Dvec = movingslope(g,2,1,data.dt);
% Dvec = Dvec - 0.04*max(Dvec); 
% %       Dvec = Dvec - mean(Dvec); 
% Dvec(Dvec < 0 )= 0;
% [up nd] = araivEst(Dvec);

% 
% up.T = zeros(size(up.trig)); 
% for i = 1:length(up.trig);
%     cut = row(up.trig(i)-40:up.trig(i));
%     DV=movingslope(cut,2,1,data.dt);
%     [m idx] = min(abs(DV));
%     up.T(i) = up.trig(i)-(40 -i);
%  
% end
% 
% 
% nd.T = zeros(size(nd.trig)); 
% for i = 1:length(nd.trig);
%     cut = row(nd.trig(i)+10:nd.trig(i)+80);
%     DV=movingslope(cut,2,1,data.dt);
%     [m idx] = min(abs(DV));
%     nd.T(i) = nd.trig(i)+10+(70 -i);
%  
% end



% valnd = zeros(size(up.T));vv=[];
% uT =[up.T n];
% for i = 1:length(uT)-1;
%  valup = row(uT(i));
%  for j=uT(i)+50:uT(i+1)-10;
%  v = row(j); 
%  vv =[vv v]; 
%  if v <= valup
%      valnd(i) = j;
%       break
%  else 
%     [m valnd(i)] = min(vv);
%       end
%  end
% end 
% 




 
 



Xnoisy = data.x/max(data.x);
denoised_norm = x/max(x);
denoised_norm = -1*denoised_norm; 

tmp = sprintf('The SNR of the denoised signal is %.2f', get_SNR(Xnoisy,denoised_norm));
disp(tmp)

% h = waitbar(0,'Plotting...');



figure(6)
subplot(4,2,1);plot(data.t,Xnoisy);
grid on
grid minor
title('Original Signal','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
xlim([min(data.t) max(data.t)]);
ax = gca;
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.05 0.75 0.420 0.150];
hold off
clear title xlabel ylabel
  

subplot(4,2,2);
[Tx, fs, as] = synsq_cwt_fw(data.t, Xnoisy, opt.nv , opt.type); 
tplot(Tx,data.t, fs); 
% spectrogram(Xnoisy,128,120,128,1E3,'yaxis'); colormap (jet)
title({'SS-CWT Spectrum'},'Rotation',0,'FontSize',13); 
xlim([min(data.t) max(data.t)]);
xlabel({'Time (s)'},'FontSize',11)
ylabel('Frequency (Hz)','FontSize',11)
ax = gca;
% ax.CLim = [-10 1000000];
ax.YAxisLocation = 'right';
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.510 0.75 0.420 0.150];
ax.CLim=[-0.00005 0.0005];
hold off
clear title xlabel ylabel


subplot(4,2,3);
plot(data.t,denoised_norm)
grid on 
grid minor
title('Denoised Signal','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
xlim([min(data.t) max(data.t)]);
ax = gca;
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.05 0.50 0.420 0.150];
hold off
clear title xlabel ylabel


subplot(4,2,4);
tplot(TT, data.t, ff); 
% spectrogram(denoised_norm,128,120,128,1E3,'yaxis'); colormap (jet)
title({'SS-CWT Spectrum'},'Rotation',0,'FontSize',13);
% xlim([min(data.t) max(data.t)]);
xlabel({'Time (s)'},'FontSize',11)
ylabel('Frequency (Hz)','FontSize',11)
ax = gca;
% ax.CLim = [-10 1000000];
ax.YAxisLocation = 'right';
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.510 0.50 0.420 0.150];
ax.CLim=[-0.35 3.5];
hold off
clear title xlabel ylabel


subplot(4,2,5);
plot(data.t,Xnoisy)
grid on 
grid minor
title('Zoomed Window- Noisy','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
xlim([min(data.t) max(data.t)]);
ax = gca;
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.05 0.28 0.200 0.150];
hold off
clear title xlabel ylabel


subplot(4,2,6);
plot(data.t,denoised_norm)
grid on 
grid minor 
title('Zoomed Window-Denoised','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
xlim([min(data.t) max(data.t)]);
ax = gca;
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.27 0.28 0.200 0.150];
hold off
clear title xlabel ylabel



subplot(4,2,7);
tplot(Tx, data.t, fs); 
% spectrogram(denoised_norm,128,120,128,1E3,'yaxis'); colormap (jet)
title({'Zoomed Spectrum'},'Rotation',0,'FontSize',13);
% xlim([min(data.t) max(data.t)]);
xlabel({'Time (s)'},'FontSize',11)
% ylabel('Frequency (Hz)','FontSize',11)
ax = gca;
% ax.CLim = [-10 1000000];
ax.YAxisLocation = 'right';
ax.YTick =[] ;
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.510 0.28 0.200 0.150];
ax.CLim=[-0.00005 0.0005];

hold off
clear title xlabel ylabel


subplot(4,2,8);
tplot(TT, data.t, ff); 
% spectrogram(denoised_norm,128,120,128,1E3,'yaxis'); colormap (jet)
title({'Zoomed Spectrum'},'Rotation',0,'FontSize',13);
% xlim([min(data.t) max(data.t)]);
xlabel({'Time (s)'},'FontSize',11)
ylabel('Frequency (Hz)','FontSize',11)
ax = gca;
% ax.CLim = [-10 1000000];
ax.YAxisLocation = 'right';
ax.TitleFontSizeMultiplier = 1.1;
ax.LabelFontSizeMultiplier=1.1;
ax.FontWeight='bold';
ax.Position=[0.73 0.28 0.200 0.150];
ax.CLim=[-0.35 3.5];
hold off
clear title xlabel ylabel





% 
% % %% Characteristic function R 
% figure (7)
% subplot 311
% plot(data.t,denoised_norm)
% grid on 
% grid minor
% title('Denoised Signal','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
% xlim([0 120]);
% 
%  ax = gca;
%  ax.TitleFontSizeMultiplier = 1.1;
%  ax.LabelFontSizeMultiplier=1.1;
%  ax.FontWeight='bold';
%  ax.Position=[0.13 0.72 0.800 0.200]; 
%  hold on
%  plot(up.T/100,0,'o','MarkerE','k','MarkerF','w')
%  
%  hold off
%  clear title xlabel ylabel
%   
% 
%  hold on 
%  subplot 312
%  plot(data.t,row);
%  grid on 
%  grid minor
%  xlim([0 120]);
%  title('Envelop Characteristic Function','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
% 
%  yrange=get(gca,'ylim');
%  for i = 1:length(out);
%  h(i) = line([up.T(i)/100,up.T(i)/100],yrange);
%  set(h(i),'Color','magenta','LineWidth',1.0);
%  hold on
%  end
% 
%  ax = gca;
%  ax.TitleFontSizeMultiplier = 1.1;
%  ax.LabelFontSizeMultiplier=1.1;
%  ax.FontWeight='bold';
%  ax.Position=[0.13 0.41 0.800 0.200];
%  hold off
%  clear title xlabel ylabel
%   
%  
% hold on 
% subplot 313
% %  plot(Dvec);
%  plot(data.t,er);
%  grid on 
%  grid minor 
%  xlim([0 120]);
%  title('Modeified Energy Ratio','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
% 
%  ax = gca;
%  ax.TitleFontSizeMultiplier = 1.1;
%  ax.LabelFontSizeMultiplier=1.1;
%  ax.FontWeight='bold';
%  ax.Position=[0.13 0.10 0.800 0.200];
%  hold on
%  plot(up.T/100,0,'o','MarkerE','k','MarkerF','w')
%  hold off
%  clear title xlabel ylabel
% 
%  
 
 
 
 
 
 
%  
%  % %% Characteristic function R 
% figure (8)
% subplot 311
% plot(data.t, denoised_norm)
% grid on 
% grid minor
% title('Denoised Signal','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
% 
% hold on 
%  yrange=get(gca,'ylim');
%  for i = 1:length(out);
%  h(i) = line([up.T(i)/100,up.T(i)/100],yrange);
%  set(h(i),'Color','magenta','LineWidth',1.0);
%  hold on
%  end
%  
%  xlim([0 120]);
%  ax = gca;
%  ax.TitleFontSizeMultiplier = 1.1;
%  ax.LabelFontSizeMultiplier=1.1;
%  ax.FontWeight='bold';
%  ax.Position=[0.13 0.72 0.250 0.200];
%  hold off
%  clear title xlabel ylabel
%   
% 
%  hold on 
% subplot 334
% 
%  plot(data.t,row);
%  grid on 
%  grid minor
%  xlim([0 120]);
%  title('Envelop Characteristic Function','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
% 
%  yrange=get(gca,'ylim');
%  for i = 1:length(out);
%  h(i) = line([up.T(i)/100,up.T(i)/100],yrange);
%  set(h(i),'Color','magenta','LineWidth',1.0);
%  hold on
%  end
%   
%  
%  ax = gca;
%  ax.TitleFontSizeMultiplier = 1.1;
%  ax.LabelFontSizeMultiplier=1.1;
%  ax.FontWeight='bold';
%  ax.Position=[0.13 0.41 0.250 0.200];
%  hold off
%  clear title xlabel ylabel
%   
% hold on 
% subplot 337
% %  plot(Dvec);
%  plot(data.t,er);
%  grid on 
%  grid minor 
%  xlim([0 120]);
%  title('Modeified Energy Ratio','Rotation',0,'FontSize',14);xlabel({'Time (s)'}); 
% 
%  ax = gca;
%  ax.TitleFontSizeMultiplier = 1.1;
%  ax.LabelFontSizeMultiplier=1.1;
%  ax.FontWeight='bold';
%  ax.Position=[0.13 0.10 0.250 0.200];
%  hold on
%  plot(up.T/100,0,'o','MarkerE','k','MarkerF','w')
%  
%  hold off
%  clear title xlabel ylabel
%  


% close(h)
toc
{'time laps' toc}

