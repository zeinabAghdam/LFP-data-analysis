clc; clear all; close all

lfpanimalDirec = '/home/sahar/Google Drive/Codes/Data_directory/Experiment_3/RatData_25April2013_6June2013/Data/Mixed-contra-ipsi/';
saveResDirec = '/home/sahar/Google Drive/Codes/GI-ICA/Application-GI-ICA/Resultsof_mainApplying_ratteData/';
%Initialization of parameters
wlen = 400;
h = 10;
nfft = 400;
subjectNumber = 1;
dataList = 2;
alpha = 1;
numSources = 10;
lastChannelNum = 12;
vLowFreq  = [1, 4, 8, 8, 12.5, 16, 32];
vHighFreq = [3.0, 7,15, 12,15.5, 31, 90 ];
assert(length(vLowFreq) == length(vHighFreq));
numBands = length(vLowFreq);
%% -------------------------------------------------
% Load the data
% -------------------------------------------------
str = funcReadFiles(  subjectNumber , dataList  );
loadDirec = strcat(lfpanimalDirec, str{1});
load(loadDirec);
trigger_signal = Y(lastChannelNum,:);
% -------------------------------------------------
% Filter LFP measurements
% -------------------------------------------------
[matLFP, samp_freq, downsampling_factor] = funcFilterRatData(Y);
% -------------------------------------------------
% Segment the LFP signals to 1.5s after the trigger
% General representation of data
% -------------------------------------------------
ht = figure;
numTrials = 60;
trialsLength = 301;
sigLength = length(matLFP);
%segmented_LFP = zeros(numSources,numTrials,trialsLength);

for nchannel = 1:numSources  
    % apply segmentation
    vecChannelSignal = matLFP(nchannel,:);
    segmented_LFP(nchannel,:,:) = ...
        funcSegmentAmplitudeData(vecChannelSignal, Y,  downsampling_factor, samp_freq);
    
    % plot different channels
    subplot(2,5,nchannel),
    plot(mean(squeeze(segmented_LFP(nchannel,:,:)))); xlim([1 length(segmented_LFP)]);
end

mystr = 'Average of evoked responses over all trials';
suptitle(mystr);
% Save the segmented LFP matrix: segmented_LFP
strName = strcat('segment_matrix_',num2str(str{1}),'_channel_',num2str(nchannel));
saveDir1 = strcat(saveResDirec,'SegmentedData/',strName);
save(saveDir1,'segmented_LFP');
% Save the plot of the segmented data
% ** ask **
%saveDir2 = strcat(saveResDirec,'fig/','plot_segmented_data_',num2str(str{1}));
%save(ht,saveDir2);
%%
% -------------------------------------------------
% Load the animal lfp data and apply bandpass filter. 
% Save the hilbert transform of the data. 
% Plot the band-passed signals of the channel.
% Save the plots.
% -------------------------------------------------
matBandPassedSignals = zeros(numSources, numBands, sigLength);
matHilbertTransform = zeros(numSources, numBands, sigLength );

for numelectrode = 1:numSources  
    h = figure;  
    for iband = 1:numBands
    
        signal = matLFP(numelectrode,:);     
        lowFreq = vLowFreq(iband);
        hiFreq = vHighFreq(iband);

        [bndpassedSignal, b] = fbandpass(signal, lowFreq, hiFreq, samp_freq);
        matBandPassedSignals(numelectrode, iband, :) = bndpassedSignal;
        matHilbertTransform(numelectrode, iband, :) = hilbert(bndpassedSignal);
        
        % plot the signal at different bands
        subplot(numBands,1,iband),
        plot(bndpassedSignal), xlim([1 length(bndpassedSignal)]);
        %hist(bndpassedSignal,30);
        titlename = strcat('band ', num2str(lowFreq),' and ', num2str(hiFreq));
        title(titlename);
        
    end
    
    titlename = strcat('ne',num2str(numelectrode));
    suptitle(titlename);
end
% save the bandpassed signals and the hilbert transformation 
saveDir3 = strcat(saveResDirec,'SegmentedData/bandpassed_subj_',num2str(str{1}));
saveDir4 = strcat(saveResDirec,'SegmentedData/hlbert_subj_',num2str(str{1}));
save(saveDir3,'matBandPassedSignals');
save(saveDir4,'matHilbertTransform');

%% -------------------------------------------------
% Apply complex ICA for every specific band on the 
%  - hilbert transformed of the signals.
% -------------------------------------------------
matReconstructedSources = zeros(iband, numSources, sigLength);
 
figure,
 for iband = 1:numBands
        iband
        S = squeeze(matHilbertTransform(:,iband,:));
        % -------------------------------------------
        % applying GI-ICA on the matrix MixedSignal.
        % data is assumed to contain noise. 
        % -------------------------------------------    
        % -------------------------------------------
        % GI-ICA : returns the reconstructed data re_dt
        %        : returns the mixingMatrix W_estimate
        %if rem(ifreq, 5) == 0
            %ifreq
        %end
        % Use random initial estimates, for now.
        A_estimate = randn(numSources)  + 1i * randn(numSources);
        A_estimate = orth(A_estimate);

        % GI-ICA
        [ re_dt, ~ ] = ComplexICA( S, A_estimate );

           % Fast CA (for debugging)
    %      meanM = repmat(mean(MixedSignal,2),1,length(MixedSignal));
    %      matCentered = MixedSignal - meanM;
    %      [E, D] = eig(cov(matCentered'));
    %      whitening_mx = sqrtm(pinv(D)) * E';
    %      whitened_dt = whitening_mx * matCentered;
    %      W = singleSourceICA(whitened_dt, 1.5);
    %      re_dt = W.' * whitened_dt;

        matReconstructedSources(iband,:,:) = re_dt;
        
        %R.matEstimatedA(ifreq,:,:) = A_estimate;
        % matIdealW(ifreq,:,:) = W_estimate;
 end

saveDir5 = strcat(saveResDirec,'SegmentedData/re_sources_',num2str(str{1}));
save(saveDir5,'matReconstructedSources');

%% Visualisation
% Plot the estimated sources at different bands.
plotBandpass(numBands, numSources, matReconstructedSources);

icounter = 1;
for iband = 1:numBands
    estimated_source = squeeze(matReconstructedSources(iband,:,:));
    for ichannel = 1:numSources
        % for every band, plot the estimate sources. 
        redt = estimated_source(ichannel,:);
        re_segmentedData = funcSegmentAmplitudeData(real(redt), Y,  downsampling_factor, samp_freq);
        
        % remove the trials with artifacts
        [re_ed_segmentedData, mat_noisyTrials] = artefact(re_segmentedData,350);
        
        subplot(numSources,numBands,icounter),
        %imagesc(re_ed_segmentedData); colormap hot;
        XX = mean(re_ed_segmentedData); plot(XX);
        xlim([1 length(XX)]);
        icounter = icounter +1;
    end
end
suptitle('Estimated Sources for numSources x numBands ')

[reconstructed_segmented] = funcSegmentAmplitudeData(matReconstructedSources, Y,  downsampling_factor, samp_freq);






