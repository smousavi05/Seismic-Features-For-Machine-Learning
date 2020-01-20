function [AICd] = akaikeC(col)
f1 = gmdistribution.fit(col,1);
f2 = gmdistribution.fit(col,2);
AICd = (f1.AIC-f2.AIC)/max([f1.AIC,f2.AIC]);