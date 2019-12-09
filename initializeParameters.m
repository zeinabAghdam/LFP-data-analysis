function [] = initializeParameters()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global windowsize
global segments_offset
global totalChannel
global totalFreq
global low
global high
global numBins
global dataDir
global subjnumber
global indStart
global indEnd

windowsize = 100;
segments_offset = 1;
totalChannel = 10;
totalFreq = 40;
low = -pi;
high = pi;
numBins = 10;
%dataDir = strcat('C:\Users\SaharAghdam\Google Drive\Codes\Data_directory\Experiment_3\RatData_25April2013_6June2013\SegmentedRatData\');
dataDir = '/home/sahar/CloudStation/RatExperiment/newSegmentedData/';
subjnumber = 1;
indStart =1;
indEnd = windowsize;

end

