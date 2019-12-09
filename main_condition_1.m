% condition 1: Filter -> Downsample -> Extract Phase -> Segment 
% In this script we segment the LFPs
% The data is bandpassed between 1Hz and 100Hz
% The LFPs are segmented into durations of 1.5 seconds 
% The segmented LFPs as well all oscillatory measurements have been saved.
clear all; close all

% Initializations
global regmode 
global nS
global lstCh
global lfpanimalDirec
global savePath
global num_ipsi 
global num_contra
global num_both 
global TFACT
global nE
global bandpass_low
global bandpass_high 
global svDirec



InitializationParams();
% The data has been separated into ipsi and contra electrodes.
% Select opt 1 for Ipis, opt 2 for contra and opt 3 for both cases. 
% L returns the number of measurements to be analyzed
opt = 3;

switch opt
    case 1
        L = num_ipsi;
    case 2
        L = num_contra;
    otherwise 
        L = num_both;
end

% In this part we set the opt to 3 as we segment all the data corresponding
% to the ipsi and contra.
for i = 1:L
    %iosition = [2,3,4,5,6,7,8,9,1,10]; % initial position of the
    %electrodes
    %vlabels = {'Infragranular e_1','e_2','e_3','e_4','e_5','e_6','e_7','Supragranular e_8','Thalamus','ALR'};
    %-----------------------------------------------------------------
    % load the data 
    % electrode position after ordering (ALR- e1- -> Thalamus)
    [data, fname] = load_data(i, opt);
    Y = data.Y;
    %Y = order_electrodes(Y,2);
    %-----------------------------------------------------------------
    
    %-----------------------------------------------------------------
    % Filter the data 
    [matLFP, samp_freq, downsampling_factor] = funcFilterRatData(Y);
    %-----------------------------------------------------------------

    %-----------------------------------------------------------------
    % change the trigger signal according to the downsampling effect 
    clear trigger_signal;
    j_out = 1;
    for j = 1:downsampling_factor:(length(Y(lstCh,:))-downsampling_factor)
        % Make sure we don't miss a trigger signal by using any
        trigger_signal(1, j_out) = any(Y(lstCh, j:(j+downsampling_factor-1)));
        j_out = j_out + 1;
    end
    %------------------------------------------------------------------
   
    %-----------------------------------------------------------------
    % Extract the phase information at different scales 
    % numChannels x numFrequencies x LengthOfSignal 
    for ichannel = 1:size(matLFP, 1)       
        signal = matLFP(ichannel,:);
        [comp_phaseData(ichannel, :, :)]...
            = extractPhase(signal,samp_freq);
    end   
    %------------------------------------------------------------------
    % The duration of data to be segmented is about 1.3 seconds, check the
    % initialization function 
    % matBT contains 
    [matSegData, matS, matSegPre, matBT] = ...
        func_SegmentData(samp_freq, matLFP, trigger_signal);
    
    % segmenting the phase matrix -> comp_phaseData
    % matSegmentedPhase ->[numFreq x nmumChannels x numSegments x Length]  
    for ifreq = 1:40 % 40: should be same as number of frequencies set in 
        pmat = squeeze(comp_phaseData(:,ifreq,:));
        
        [matSegmentedPhase(ifreq,:,:,:), matP(ifreq,:,:),~,~] = ...
            func_SegmentData(samp_freq, pmat, trigger_signal);
        
    end
    %------------------------------------------------------------------
    
    %------------------------------------------------------------------
    % For every channel we remove the trials with a trial above 300 
    % mat_ar_sweeps stands for matrix of artefact removed sweeps
    for ich = 1:nE
       temp = squeeze(matSegData(ich,:,:));
       [mat_tmp_sweeps , ~, nsweeps(ich)] = artefact(temp);
       if ~isempty(nsweeps(ich))
           if nsweeps(ich) < min(nsweeps(1:ich-1))
               temp2 = mat_ar_sweeps(:,1:nsweeps(ich),:);
               clear mat_ar_sweeps
               mat_ar_sweeps(:,:,:) = temp2;
               mat_ar_sweeps(ich,:,:) = mat_tmp_sweeps(1:nsweeps(ich),:);
           else
               minLen = min(nsweeps);
               mat_ar_sweeps(ich,:,:) = mat_tmp_sweeps(1:minLen,:);
           end
       end
           
       clear mat_tmp_sweeps
       clear temp 
    end
    % ---------------------------------------------------------------
    % In this part, we compute the phase after segmentation 
    % The phase is computed for every trial 
    % phaseMatPerSeg -> numChannels x numFreq x numTrials x T
    % ----------------------------------------------------------------
    [phaseMatPerSeg, phaseMatPerSegComplex]  = compute_phase_PerSegment(mat_ar_sweeps, samp_freq);
    saveResults(matSegData, matS, mat_ar_sweeps, matSegmentedPhase, phaseMatPerSeg, phaseMatPerSegComplex,matSegPre, matBT, fname)
    %------------------------------------------------------------------  
    clear matSegData
    clear matS
    clear mat_ar_sweeps
    clear noisyTrials
    clear nsweeps 
    clear phaseMatPerSeg
    clear matSegmentedPhase
    clear matP
    
end

