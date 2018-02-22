function [adjMat] = getAdjMat(pos,plotGraph,roiVals)
    if nargin < 2
        plotGraph = 0;
    end
    N = size(pos,1);
    adjMat = zeros(N,N);
    for ii = 1:N
       for jj = 1:N
            adjMat(ii,jj) = norm(pos(ii,:) - pos(jj,:)); 
       end
    end
    if plotGraph
        tmpMat = cov(normc(roiVals)); 
        g = graph(tmpMat.*(tmpMat > 1.2*10^-6),'OmitSelfLoops'); 
        LWidths = 5*abs(g.Edges.Weight)/max(abs(g.Edges.Weight)); 
        plot(g ,'-.dr', 'YData' , -pos(:,1) ,...
            'XData' , pos(:,2),'EdgeLabel',...
            g.Edges.Weight,'LineWidth',LWidths)
    end
end