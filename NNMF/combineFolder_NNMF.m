function [] = combineFolder_NNMF(pathTo)
    close all
    tmp = ls(pathTo);
    firstExperiment = 1;
    fileList = strsplit(tmp);
    fileList = fileList(1:end-1);
    allMeans = [];
    figure();
    colormap(summer); shading flat;
    for ii = 1:length(fileList)

        subplot(4,4,2*ii-1)
        file = fileList{ii}
        [~ , means , ~ ,  ~ , ~] = importCSV(fullfile(pathTo,file));
        [W,H] = nnmf(transp(means),3);
        plot(transp(H))
        xlim([0,size(H,2)])
        if ii == 1
        legend({'Factor 1','Factor 2','Factor 3'});%,'Factor 4',...
        %'Factor 5','Factor 6','Factor 7','Factor 8'})
        end
        subplot(4,4,2*ii)
        imagesc(W);colorbar;
        allMeans = [allMeans ; transp(means)];
    end
    figure()
    subplot(2,1,1)
    [W,H] = nnmf(allMeans,8);
    plot(transp(H));        xlim([0,size(H,2)])
    legend({'Factor 1','Factor 2','Factor 3','Factor 4',...
        'Factor 5','Factor 6','Factor 7','Factor 8'})
    xlabel('Time t'); ylabel('dF/F0');
    subplot(2,1,2)
    imagesc(W);colorbar;
end
