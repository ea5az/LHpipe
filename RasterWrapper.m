%%
% A wrapper for the plotSpikeRaster tool from the internet. If provided
% with participation rates, also plots network events.
%   Input: binMat - a binary matrix of spike events
%          parRates - a vector of participation rates per time step
%   Output: f - figure handle to scatter plot
%
function [f] = RasterWrapper(binMat,evStarTimes,evDurations,evRates,ParRate)
    addpath(genpath('./tools/plotRaster'));
    %f = figure();
    hold on;
    % If participation rates available
    if nargin > 1
       yyaxis right
       % Only participation > 20% is network event, Siegel 2012
       scatter(evStarTimes , ones(length(evStarTimes),1),'*')
       scatter(evStarTimes+evDurations' , ones(length(evStarTimes),1),'*')

       for ii = 1:length(evRates)
           text(evStarTimes(ii) , 0.9 , sprintf('%1.2f',evRates(ii)),'Rotation',290)
           plot([evStarTimes(ii) evStarTimes(ii)+evDurations(ii)] , [1 1] , '-')
       end
       xlim([0 , length(binMat)])
       set(gca , 'ytick' , [])
       yyaxis left
    end
    plotSpikeRaster(logical(transp(binMat(:,sum(binMat,1)>0))),'PlotType','imagesc');
    hold on
    plot(ParRate*size(binMat,2))
end