function [ segment_matrix ] = funcAmplitudeSegmentation( Y, vecChannelSignal,downsampling_factor )
global samp_freq
% Segementation of data according to trigger / generation of data cube
disp('Stage 5 - Data segmentation');
disp('       Trigger downsampling');

clear trigger_signal;
i_out = 1;
for i = 1:downsampling_factor:(length(Y(12,:))-downsampling_factor)
    % Make sure we don't miss a trigger signal by using any
    trigger_signal(1, i_out) = any(Y(12, i:(i+downsampling_factor-1)));
    i_out = i_out + 1;
end

trigger_count = 0;
for i4 = 1:length(trigger_signal)-1;
    if trigger_signal(1,i4) == 0 && trigger_signal(1,i4+1) == 1
        trigger_count = trigger_count + 1;
    end
end

%trigger_count
disp('   ---       Data segmentation');
trial_dur = 1.5; % 1.5 in seconds
num_samples_per_trial = trial_dur * samp_freq;
clear segment_matrix;
row_count = 0;
for i5 = 1:length(trigger_signal)-num_samples_per_trial
    if trigger_signal(1,i5) == 0 && trigger_signal(1,i5+1) == 1
        row_count = row_count + 1;
        segment_matrix(row_count,:) = vecChannelSignal(i5:i5+num_samples_per_trial);
    end
end

end

