function [mat_segmentedData,mat_S, mat_bfT_segmentedData,mat_bT ] =...
    func_SegmentData(samp_freq, matLFP, trigger_signal)

global trial_dur

num_samples_bf_trial = round(0.1 * samp_freq); % 30ms before the stimulus
num_samples_per_trial = round(trial_dur * samp_freq);

for iChannel = 1:size(matLFP,1)
    row_count = 0;
    bTrigger = false;
    vecChannelSignal = matLFP(iChannel,:);

    for ix = 1:length(trigger_signal)-num_samples_per_trial

        if trigger_signal(1,ix) == 0 && trigger_signal(1,ix+1) == 1
            tmp = ix+num_samples_per_trial;
            if tmp <= length(vecChannelSignal)
                row_count = row_count + 1;
                segment_matrix(row_count,:) = vecChannelSignal(ix:ix+num_samples_per_trial);
                if bTrigger
                    % responses before the trigger 
                    segment_bt_matrix(row_count,:) = vecChannelSignal(ix-num_samples_bf_trial:ix);
                end
            end
            bTrigger = true; 
        end
    end
    mat_segmentedData(iChannel,:,:) = segment_matrix;
    mat_S(iChannel,:) = reshape(segment_matrix, 1, size(segment_matrix,1)*size(segment_matrix,2));
    
    % segment before trigger data 
    mat_bfT_segmentedData(iChannel,:,:) = segment_bt_matrix;
    mat_bT(iChannel,:) = reshape(segment_bt_matrix, 1, size(segment_bt_matrix,1)*size(segment_bt_matrix,2));
    
end


end