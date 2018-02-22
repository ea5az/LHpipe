%%%
% Extracts participation rate from a binary matrix of spike events and a
% provided size of a window, default = 3 (Siegel 2012).
%   Input: binMat - a binary matrix of spike events
%          windowSize - size of the window, this delay is considered
%          simultaneous
%   Output : ParRate : A vector with participation rates per time step
%

function [ParRate , duds] = get_ParRate(binMat , windowSize , dudslim)
    if nargin < 3
        dudslim = 20;
    end
    if nargin  < 2
        windowSize = 3;
    end
    M = size(binMat,1);
    N = size(binMat,2);
    ParRate = [zeros(1,(windowSize-1)/2)];
    % Exclude ROIs that don't produce any events from the analysis. The
    % marker probably didn't reach these.
    duds = sum(sum(binMat,1)<dudslim);

    %disp(sprintf('Dudsrate %f',duds/N))
    
    for ii = 1+(windowSize-1)/2:M-(windowSize-1)/2
        window = transp(binMat(ii-(windowSize-1)/2:ii+(windowSize-1)/2,:));
        parWindow = sum(window,2) > 0 ;
        ParRate = [ParRate , sum(parWindow)/(N-duds)];
    end
    ParRate = [ParRate zeros(1,(windowSize-1)/2)];
end