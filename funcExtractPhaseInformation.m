function [w_phase] = funcExtractPhaseInformation(w_scales,y)

w_coefs = [];
w_coefs(:,:) = cwt(y,w_scales,'cgau6');  

%%%%% Extraction of instantaneous phase %%%%%
disp('Stage 4 - Phase extraction');
for i = 1:length(w_scales)
    w_phase(i,:) = angle(w_coefs(i,:));
end;
%%%%% Save Data %%%%%
disp('Saving phase Data ....  ')
%save_nam = (['phase_angles_scales_3-1000_20_12_E' num2str(chan_select)]); % Adapt file name

end
