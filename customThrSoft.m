function [ dnCuCoef ] = customThrSoft(wlCoef,opt,nn)

% Estimate the unuversal thresholding level.


for k =1:nn.ny;
Wx_fine = abs(wlCoef(k, 1:nn.nu(1)));
gamma = sqrt(2*log(nn.n)) * mad( abs(Wx_fine (:))) * 1.4826;
%
% gamma = thselect(wlCoef,'sqtwolog')

% Customized Thresholding
lam = 0.9*gamma; % cutoff value
al = 0;         % shape parameter

% for j = 1:nn.na
for i=1:nn.n
    val = abs(wlCoef(k,i));
   if val <= lam
      wlCoef(k,i) = 0;

   elseif val >= gamma
      wlCoef(k,i) = wlCoef(k,i)-sign(wlCoef(k,i))*(1-al)*gamma;

   else
      wlCoef(k,i) = al*gamma*((val-lam)./(gamma-lam))^2*(al-3)*((val-lam)./(gamma-lam))+4-al;
     end
end
end
dnCuCoef = wlCoef;





