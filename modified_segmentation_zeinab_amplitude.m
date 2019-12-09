%% Extraction of scale and instantaneous
%% phase before segmentation.
%clc; clear all; close all;
%disp('Loading data ....');
%load('rn_20122012_03_ipsi'); % Insert filename here
%% 1Hz to 300 Hz bandpass filtering and removal of 50Hz and simple
%% harmonics
clear all
%%%%% Variable declaration %%%%%
global samp_freq
no_of_channels = 12;
filtered_Y2 = [];%zeros(no_of_channels,length(Y));
filtered_Y3 = [];%zeros(no_of_channels,length(Y));
filtered_Y4 = [];
segment_matrix = [];

% load data 
subject_number = 1;
str = funcReadFiles( subject_number , 1  );
load(strcat('C:\Users\SaharAghdam\Google Drive\Codes\Data_directory\Experiment_3\RatData_25April2013_6June2013\Data\',str{1}));

%samp_freq = 1920; %9600;
samp_freq = length(Y)/Y(1,end);
%samp_freq = 4800;

filter_order = 2000;
no_of_iterations = 5;

%%%%% Bandpass filtering and downsampling %%%%%
disp('Stage 1 - Bandpass filtering and downsampling');
nyquist_freq = samp_freq/2;
bandpass_low = 2; % Herz
bandpass_high = 300; % Herz
cut_off = [(bandpass_low/nyquist_freq) (bandpass_high/nyquist_freq)];
b = fir1(filter_order,cut_off,'bandpass');

%%
% conv with a Gaussian filter 
% Construct blurring window.
% After convolving the signals with a Gaussian filter, apply the filtering.
% windowWidth = int16(5);
% halfWidth = windowWidth / 2;
% gaussFilter = gausswin(5);
% gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
% 
% for count_a2 = 2:no_of_channels-1
%     % Do the blur.
%     smoothed_Y(count_a2-1,:) = conv(Y(count_a2,:),gaussFilter); 
%     Y_filt(count_a2-1,:) = detrend(smoothed_Y(count_a2-1,:),'constant');
%     filtered_Y2(count_a2-1,:) = filtfilt(b,1,Y_filt(count_a2-1,:));
% end
%%
% Uncomment here if the section above is not used.
for count_a = 2:no_of_channels-1;
    Y(count_a,:) = detrend(Y(count_a,:),'constant');
    filtered_Y2(count_a-1,:) = filtfilt(b,1,Y(count_a,:));
end;

%%
disp('Downsapling ....')
% Downsample to twice the new Nyquist frequency
downsampling_factor = floor(samp_freq/(bandpass_high * 2));
%downsampling_factor = 1; % Uncommend to disable downsampling
% The channels are from 2 to 11
for count_b = 1:size(filtered_Y2,1);
    filtered_Y4(count_b,:) = downsample(filtered_Y2(count_b,:),downsampling_factor);
end
samp_freq = samp_freq / downsampling_factor;
nyquist_freq = samp_freq / 2;
%%
%%%%% Removal of 50Hz (2Hz bandwidth) and simple harmonics %%%%%
disp('Stage 2 - Filtering 50Hz and harmonics');
filter_order_stop = 2000;

for filter_count = 1:no_of_iterations
    cut_off_stop = [(filter_count*50-1)/nyquist_freq (filter_count*50+1)/nyquist_freq];
    b2 = fir1(filter_order_stop,cut_off_stop,'stop');

    for count_a2 = 1:size(filtered_Y4,1);
        filtered_Y4(count_a2,:) = filtfilt(b2,1,filtered_Y4(count_a2,:));
    end;
end;
%%%%% Save data %%%%%
disp('Saving data');
%save('filtered_raw_data_all_channels_rat_20_12.mat','filtered_Y4');
%%
for i_chan_select=1:size(filtered_Y4,1)
    segment_matrix = funcAmplitudeSegmentation(Y,filtered_Y4(i_chan_select,:),downsampling_factor);
    subplot(5,2,i_chan_select),
    imagesc(segment_matrix), 
end
