function [denoised ] = synchDenois(dec,opt,data,nn);


% soft thresholding long periods
[na, n] = size(dec.long.wl);
gamma = sqrt(2*log(n)) * mad( abs(dec.long.wl(:))) * 1.4826;
dec.long.wl= SoftThresh(dec.long.wl,gamma);

% Synchrosqueezing
[dn.long ,fs.long] =sst(dec.long,opt);
[Tx, ff, as] = synsq_cwt_fw(data.t, data.x, opt.nv , opt.type); 

[tx.rest ,fs.rest] =sst(dec.rest,opt);
[tx.noise ,fs.noise] =sst(dec.noise,opt);

gammaN = sqrt(2*log(n)) * mad( abs(tx.noise(:))) * 1.4826;
highN = HardThresh(tx.noise,gammaN);
sigma =  mad(highN(:))./0.6745;
Mmax = mean(max(abs(highN)));
Sig = mad(tx.noise(:))./0.6745;


% finding the frequency band of highest noise concentration 
[nu.noise nd.noise]= sep(highN ,fs.noise, dec.noise.t );

fUp = []; fDn = [];
for i = 1:length(nu.noise);
   fUp = [ fUp fs.noise(nu.noise(i)) ];
   fDn = [ fDn fs.noise(nd.noise(i)) ];
end


idUpRest = []; idDnRest = [];
for j = 1:length(fUp);
tmp = abs(fs.rest(:) - fUp(j));
[m idx] = min(tmp) ;
idUpRest = [idUpRest idx];

tmp = abs(fs.rest(:) - fDn(j));
[m idx] = min(tmp) ;
idDnRest = [idDnRest idx];

end 


% taking care of the noise
[na n] = size(tx.rest);
for i = 1:length(idUpRest);

for j=idUpRest(i):idDnRest(i);
for k = 1:n;
    val = abs(tx.rest(j,k));
    
  if abs(val) > Mmax  
    
    res = (abs(val) - Mmax);
    tx.rest(j,k) = tx.rest(j,k)*(res/val);
    
  elseif abs(val) <= Mmax &  abs(val) > gammaN;
      
    tx.rest(j,k) = 0;
  else 
     tx.rest(j,k) = tx.rest(j,k)/gammaN;
  end  

end
end 
end


% asembelyy 
dnFnl = zeros(nn.na, nn.n);
dnFnl(nn.ny+1:nn.na,:)= dn.long;
dnFnl(1:nn.ny,1:nn.n)= tx.rest;

dnFnl(isnan(dnFnl)) = 0;
denoised = synsq_cwt_iw(dnFnl, ff, opt);



% figure(79);plot(denoised)

% for k = 1 : length(idUpRest)
%  tx.rest(idUpRest(k):idDnRest(k),:) = 0;
% end