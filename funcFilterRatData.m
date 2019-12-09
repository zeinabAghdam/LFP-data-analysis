% Filtering the LFP measurements between 1Hz and 100Hz, removing the 50Hz
% harmonics 
function [filtered_Y3, samp_freq, downsampling_factor] = funcFilterRatData(Y)  

global no_of_channels 
global bandpass_low
global bandpass_high 

Y1 = zeros(no_of_channels-2,length(Y));
filtered_Y = zeros(no_of_channels-2,length(Y));
filtered_Y2 = zeros(no_of_channels-2,length(Y));

samp_freq = length(Y)/Y(1,end);
%filter_order = 500; 
filter_order = 1000;
nyquist_freq = samp_freq/2;

%bandpass_low = 2; % Hertz
%bandpass_high = 100; % Hertz

cut_off = [(bandpass_low/nyquist_freq) (bandpass_high/nyquist_freq)];
b = fir1(filter_order,cut_off,'bandpass');

for count_a = 2:no_of_channels-1
    Y1(count_a-1,:) = detrend(Y(count_a,:),'constant');
    filtered_Y(count_a-1,:) = filtfilt(b,1,Y1(count_a-1,:));
end

clear b
[b,a] = iircomb(round(samp_freq/50),50/nyquist_freq/3,'notch');

%fs = samp_freq; fo = 50;  q = 35; bw = (fo/(fs/2))/q;
%[b,a] = iircomb(floor(fs/fo),bw,'notch'); % Note type flag 'notch'
%fvtool(b,a);

% removing 50Hz 
disp('removing 50 Hz')
for count_a2 = 1:size(filtered_Y,1)
   filtered_Y2(count_a2,:) = filtfilt(b,a,filtered_Y(count_a2,:));
end

filtered_Y2 = filtered_Y2(:,2*filter_order:end-2*filter_order);

% Downsampling
disp('Downsampling data ... ')
downsampling_factor = floor(samp_freq/(bandpass_high * 2));
filtered_Y3 = zeros(no_of_channels-2,ceil(length(filtered_Y2)/downsampling_factor));
for count_b = 1:size(filtered_Y2,1);
    filtered_Y3(count_b,:) = interp1(filtered_Y2(count_b,:),1:downsampling_factor:length(filtered_Y2),'spline');
    %filtered_Y3(count_b,:) = downsample(filtered_Y2(count_b,:),downsampling_factor);
end

samp_freq = samp_freq / downsampling_factor;
%nyquist_freq = samp_freq / 2;
end
