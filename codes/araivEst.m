function [ up nd ] = araivEst(f)

% finding beginig of triggers
n = length(f);
out = [];

for i = 2:n-1
if (f(i-1) == 0 & f(i+1)> 0 & f(i) == 0);
    out = [out i+1];
    
  end
end
up.all = out; 


pick = [out(1)];
for j = 2:length(out);
     dis = out(j) - out(j-1);
  if  dis > 100;
    pick = [pick out(j)];
  end 
end
% 
% pick = pick - 100;
up.trig = pick;




% now end of trigger
out = [];
for i = 2:n-1
if (f(i-1) > 0 & f(i+1) == 0 & f(i) == 0);
    out = [out i+1];
    
  end
end

nd.all = out;



pick = [out(1)];
for j = 2:length(out);
     dis = out(j) - out(j-1);
  if  (abs(dis)>100);
    pick = [pick out(j)];
  end 
end

% pick = pick + 200;
nd.trig = pick;


