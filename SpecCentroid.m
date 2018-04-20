function [C, CMean, CSD, CMax] = SpecCentroid(snd,fs,nfft,window,noverlap)

%CENTROID Calculates spectral centroid of a sound.
%	C = CENRTROID(SOUND, FS) calculates the spectral centroid of a given
%	sound. Fs is the sample frequency; the window length is the minimum of 2048 and
%	the sample length; overlap is 80%.
%   
%   [C, CMean, CSD, CMax] = CENRTROID(SOUND, FS) also calculates the mean,
%   standard deviation and maximum of spectral centroid.
% 
%   C = CENRTROID(SOUND,FS,NFFT); 
%   C = CENRTROID(SOUND,FS,NFFT,WINDOW);
%   C = CENRTROID(SOUND,FS,NFFT,WINDOW,NOVERLAP)
%   Read the manual of the SPECGRAM command for the parameters.
%
% Example:
%	[s, fs] = wavread('testsound.wav');
%	c = centroid(s, fs);
% Uses:
%	Matlab Signal Processing Toolbox
%
% References:
%	Tzanetakis, G., Essl, G. and Cook, P.
%	In: Proc. Int. Symposium on Music Information 
%	Retrieval (ISMIR), Bloomington, Indiana, 2001 
%
% Frederik Nagel and Michael Groﬂbach
% Institute of Music Physiology and Musicians' Medicine
% Hanover University of Music and Drama 
% Hannover
% Germany
%
% e-mail: frederik.nagel@hmt-hannover.de
% homepage: http://www.immm.hmt-hannover.de
%
% May 29, 2006.
%
% See also SPECGRAM, FFT

error(nargchk(2, 5, nargin))
if(nargin==2)
    nfft = min([length(snd) 2048]);
    window = nfft;
    noverlap = round(window*.8);
    s = specgram(snd, nfft, fs, window, noverlap);
elseif (nargin==3)
    s = specgram(snd, nfft, fs);
elseif (nargin==4)
    s = specgram(snd, nfft, fs, window);
elseif (nargin==5)
    s = specgram(snd, nfft, fs, window, noverlap);    
end

C = sum((repmat((1:size(s,1))',1,size(s,2)) .* abs(s))) ./ sum(abs(s));
CMean = mean(C);
CSD = std(C);
CMax = max(C);