function E = feature_energy(window)

% function E = feature_energy(window);
%
% This function calculates the energy of an audio frame.
% 
% ARGUMENTS:
% - window: 	an array that contains the audio samples of the input frame
% 
% RETURN:
% - E:		the computed energy value
%
% (c) 2014 T. Giannakopoulos, A. Pikrakis

E = (1/(length(window))) * sum(abs(window.^2));
