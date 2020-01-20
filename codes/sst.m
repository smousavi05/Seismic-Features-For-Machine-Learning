function [Tx ,fs] =sst(wev,opt)

%% Synchrosqueezing 
% Approximate instantaneous frequency
w = imag(wev.dWx ./ wev.wl / (2*pi));
 
% Calculate the synchrosqueezed frequency decomposition
[Tx,fs] = synsq_cwt_squeeze(wev.wl, w, wev.t, opt.nv, opt);