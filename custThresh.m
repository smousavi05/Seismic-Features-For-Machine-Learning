function [ denoised opp] = custThresh(wlCoef,opt,data)


% Estimate the unuversal thresholding level.
[na, n] = size(wlCoef);
Wx_fine = abs(wlCoef(1:opt.nv,1:n));
gamma = sqrt(2*log(n)) * mad( abs(Wx_fine (:))) * 1.4826;



x=[]; y=[];num=0;op=[];
for k =0.1:0.1:0.9
    
% Customized Thresholding 
lam =  k*gamma; % cutoff value

for m =0:0.2:1;
al = m;         % shape parameter
num =num+1; 

wdn = wlCoef;
for j = 1:na
for i=1:n
    val = abs(wdn(j,i));
   if val <= lam
      wdn(j,i) = 0;
        
   elseif val >= gamma
      wdn(j,i) = wdn(j,i)-sign(wdn(j,i))*(1-al)*gamma;
        
   else
      wdn(j,i) = al*gamma*((val-lam)./(gamma-lam))^2*(al-3)*((val-lam)./(gamma-lam))+4-al;    
     end
end
end
wdn(isnan(wdn)) = 0;
xx = cwt_iw(wdn, opt.type, opt, opt.nv);
x{num}=xx;
op{num}=[k m];
RMSE = sqrt(mean((xx' - data.x).^2));
y=[y RMSE];

end
end

[n idx]=min(y);
denoised = x{idx};
opp = op{idx};
