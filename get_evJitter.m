%%
% extracts the jitter from a given array of events. Jitter is defined as
% the standard deviation of spiking times within the event. If a neuron
% fired several times during the event, the average over these spikes is
% taken.
%   Input:   raster - binary matrix of spike events
%            evStartTime - array of starting times
%            evDurations - array of event durations
%   Output:  evJitter - array of jitter per event
%
function evJitter = get_evJitter(raster , evStartTime , evDurations)
    N = size(raster,2); 
    evJitter = zeros(length(evStartTime),1);
    for ii = 1:length(evStartTime)
        currStart = evStartTime(ii);
        if evDurations(ii) + currStart > size(raster,1)
            window = raster(currStart:end,:);
        else
            window = raster(currStart:currStart+evDurations(ii),:);
        end
        M = size(window,1);
        % convert binary matrix to matrix where nonzero indices indicate
        % spike time of corresponding event
        window = repmat(transp(1:M),1,N).*window;
        window(window == 0) = NaN;
        eventMean = nanmean(window(:));
        neuroMean = window(~isnan(window));%nanmean(window,1);
        %neuroMean = neuroMean(~isnan(neuroMean));
        % standard deviation
        evJitter(ii) = sqrt( sum((neuroMean - eventMean).^2)/(length(neuroMean)-1)  );
    end
end
