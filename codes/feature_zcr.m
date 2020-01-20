function Z = feature_zcr(window);

% function  Z = feature_zcr(window);
%
% This function calculates the zero crossing rate of an audio frame.
% 
% ARGUMENTS:
% - window: 	an array that contains the audio samples of the input frame
% 
% RETURN:
% - Z:		the computed zero crossing rate value
%

window2 = zeros(size(window));
window2(2:end) = window(1:end-1);
Z = (1/(2*length(window))) * sum(abs(sign(window)-sign(window2)));



