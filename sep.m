function   [nu nd]= sep(w,fs,t)


[na n]=size(w);col=0;
for i = 1:n
    col = col+abs(w(:,i));
end
xsmooth = smooth(col,0.05,'loess');

msmooth = xsmooth - mean(xsmooth);
msmooth(msmooth<0) = 0;

nu = []; nd = [];
for i =2:length(msmooth)-1;
    if (msmooth(i-1) == 0 & msmooth(i+1) > 0 & msmooth(i) == 0);
        nu = [nu i+1];
    end
    if (msmooth(i-1) > 0 & msmooth(i+1) == 0 & msmooth(i) == 0);
        nd = [nd i+1];
    end
end


% 
% figure
% 
%  subplot 121
%  hold on
%  a = 1:length(col);
%  a = a';
%  bar(a,col);
%  xlim([0 length(col)]);
%  title({'Stacked';'Coefficients'});
%  ylabel({'Magnitude'});
% 
%  hold on
%   yrange=get(gca,'ylim');
%   
%  for i = 1:length(nu);
%  
%  h = line([nu(i),nu(i)],yrange);
%  set(h,'Color','magenta','LineWidth',1.5);
%  
%   h = line([nd(i),nd(i)],yrange);
%  set(h,'Color','magenta','LineWidth',1.5);
%  
%  end
%  hold off
%   
% ax = gca;
% ax.XTick = [];
% ax.TitleFontSizeMultiplier = 1.1;
% ax.LabelFontSizeMultiplier=1.1;
% ax.FontWeight='bold';
% ax.Position=[0.15 0.40 0.100 0.300];
% ax.View=[270 90]
% hold off
% clear title xlabel ylabel
% 
% 
% subplot 122
% hold on
% tplot(w, t, fs);
% title('SS-CWT of Noise ');
% xlabel('Time (S)');
% ylabel('Frequency (Hz)');
% 
% ax = gca;
% ax.Position=[0.265 0.40 0.600 0.300];
% ax.YAxisLocation = 'right';
% ax.CLim=[-10 100];
% ax.FontSize=12;
% ax.TitleFontSizeMultiplier = 1.3;
% ax.LabelFontSizeMultiplier=1.1;
% ax.FontWeight='bold';
% hold off 












 
 
 
 

 
 
 
 

