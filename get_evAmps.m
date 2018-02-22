%%
% extracts event amplitudes from given recorded values from ROI by
% averaging over the maximal amplitude during the duration of the event
%   Input: roiVals - time series of measured dF/F0 signal
%          evStartTime - start time of events
%          evDurations - durations of events
%   Output: evAmps - Amplitudes of events
%   
function evAmps = get_evAmps(roiVals  , raster ,  evStartTime , evDurations)
    evAmps = zeros(length(evStartTime),1);
    for ii = 1:length(evStartTime)
        try
            window = roiVals(evStartTime(ii):evStartTime(ii)+ evDurations(ii),:);
            rWindow = raster(evStartTime(ii):evStartTime(ii)+ evDurations(ii),:);
        catch
            window = roiVals(evStartTime(ii):end,:);
            rwindow = raster(evStartTime(ii):end,:);
        end
        mask = sum(rWindow,1) > 0;
        rWindow = window(:,mask);
        evAmps(ii) = mean(max(rWindow));
    end
end