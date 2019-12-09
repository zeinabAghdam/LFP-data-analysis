function [w_phase] = extractPhase(y,Fs)
    % the settings here should be the same as in extractPhase
    
global bandpass_low
global bandpass_high 

% compute the phase information at 40 different scales. 
% maybe set better values later. (TODO)
% the frequency should be between bandpass_low and bandpass_high
% this will cover between 1.9 (~2Hz) - 97.7 (~98Hz) 
% TODO: new version of cwt different!!!
w_scales = logspace(0.01, 1.72 ,  40); 
% these are the range of frequencies it covers
w_freq = scal2frq(w_scales, 'cgau6', 1/Fs); 

% compute the wavelet coefficients for all the scales 
w_coefs = [];
w_coefs(:,:) = cwt(y,w_scales,'cgau6');  

for i = 1:length(w_freq)
    w_phase(i,:) = angle(w_coefs(i,:));
end

            
end
