function xf=hipass(x,Fs,co,npol,npas,tipe,trending)
% xf=HIPASS(x,Fs,co,npol,npass,tipe,trending)
%
% Filters signal 'x' with filter 'tipe' and corner
% frequency 'co' in Hz with 'npol' number of poles and in 
% 'npas' passes. Sampling rate is 'Fs' in Hz.
%
% INPUT:
%
% x         The signal
% Fs        Its sampling frequency (default: 110)
% co        The corner frequency, in Hz (default: 5)
% npol      The number of poles (default: 2)
% npas      The number of passes (default: 1)
% tipe      The filter name (default: 'butter')
% trending  'linear' or 'constant' (default: 'linear')
%
% Compare in SAC hp butter co 5 n 2 p 1
% See FilterTest for comparison.
%
% See also LOWPASS, BANDPASS
%
% Last modified by fjsimons-at-alum.mit.edu, 09/19/2007

defval('npol',2)
defval('npas',1)
defval('co',5)
defval('Fs',110)
defval('tipe','butter')
defval('trending','linear')

disp(sprintf('HIPASS %3.3f Hz %i pass %i poles %s',co,npas,npol,tipe))
						
% Corner frequency is in Hertz, now it is as a fraction of
% half the sampling rate.
Wn=2*co/Fs;

[B,A]=feval(tipe,npol,Wn,'high');

xf=filter(B,A,detrend(x(:),trending));

if npas==2
  xf=flipud(filter(B,A,detrend(flipud(xf(:)),trending)));  
end
