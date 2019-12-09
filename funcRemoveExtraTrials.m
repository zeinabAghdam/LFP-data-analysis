function [myData] = funcRemoveExtraTrials(tmp,vecSweepIndices)

% return the new matrix of phase which does not contain the trials in which
% the amplitude data was above the threshold. 
total_freqNumber = size(tmp,2);

for ix = 1:total_freqNumber
    matData = squeeze(tmp(:,ix,:));
    matData(vecSweepIndices,:) = [];
    myData(:,ix,:) = matData;
end

end

