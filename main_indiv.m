tab = savTab;
params = getParams(DATASET); flags = getFlags(DATASET);

pathTo = params.TRACESpath;
fileList = getFileList(pathTo);
fixThresh = params.fixThresh;
% Siegel proposes to cutoff at 20%
LOWERBOUND = params.LOWERBOUND;
lAmps = [];
hAmps = [];

cID = NaN;
for ii = 1:length(fileList)
    fileStr = fileList{ii}
    % extract expriment identifiers
    [expID , expID2] = getExpID(fileStr,params);

    % get the average activity from the ROI's.
    [readIn ,skip] = importCSV([pathTo fileStr] , params , flags);
    if skip
        continue
    end
    if cID ~= expID || logical((cCond ~= readIn.cond)*flags.changeAtCondChange)
        [thrTraces , thrDiffTraces] = get_Threshold(readIn.means,fixThresh,flags); 
        cID = expID;
        cCond = readIn.cond;
    elseif logical((cCond ~= readIn.cond)*~flags.changeAtCondChange)
        sc1 = mean(thrTraces); sc2 = mean(thrDiffTraces);
        [thrTraces , thrDiffTraces] = get_Threshold(readIn.means,fixThresh,flags);
        thrTraces = thrTraces*sc1/mean(thrTraces); thrDiffTraces = thrDiffTraces*sc2/mean(thrDiffTraces);
        cID = expID;
        cCond = readIn.cond;
    end
    
    [raster , ~] = TracesToSpikeTimes(readIn.means , thrTraces , params.trDiffRat*thrDiffTraces);

    sTab = tab(tab.id == expID & tab.id2 == expID2,:);
    lAcc = []; nArr = [];
    for jj = 1:size(readIn.means , 2)
        for kk = 1:length(sTab.startTimes)
            cStartTime = sTab.startTimes(kk);
            cDur = sTab.durations(kk);
            sRaster = raster(cStartTime:cStartTime+cDur , jj);
            sMean = readIn.means(cStartTime:cStartTime+cDur , jj);
            if sum(sRaster) > 0
                if sTab.rates(kk) > params.lEventUpper
                    lAmps = [lAmps nanmean(lAcc)]; hAmps = [hAmps nanmax(sMean)];%/nanmean(nArr)
                    lAcc = [];
                else
                    if sTab.rates(kk) > 0.4 & sTab.rates(kk) < 0.7
                        nArr = [nArr nanmax(sMean)];
                    end
                    lAcc = [lAcc nanmax(sMean)];
                end
                
            end
        end
    end
end
%%
rhAmps = hAmps(~isnan(lAmps));
rlAmps = lAmps(~isnan(lAmps));
[R,P,RL,RU] = corrcoef(rlAmps,rhAmps);
ccoeff = R(1,2);
figure()
scatter(rlAmps,rhAmps ,45,'MarkerFaceColor',rgb('gray'),'MarkerEdgeColor',rgb('black'))
lsline
%colorbar
title(sprintf('correlation: %1.2f',ccoeff))
xlabel('Mean L-event amplitude ^{F}/_{F0}')
ylabel('H-event amplitude ^{F}/_{F0}')

function fileList = getFileList(pathTo)
    tmp = ls(pathTo);
    fileList = strsplit(tmp);
    fileList = fileList(1:end-1);
    fileList = sort(fileList);
     normidx = find(strcmp(fileList,'norm'));
    if ~isempty(normidx) fileList(normidx) = []; end
    matidx =  find(strcmp(fileList,'mat'));
    if ~isempty(matidx) fileList(matidx) = []; end
    arxidx =  find(strcmp(fileList,'arx'));
    if ~isempty(arxidx) fileList(arxidx) = []; end
end


function  [expID , expID2] = getExpID(fileStr , params)
        % add experiment ID from string
        %expIDstr = file{1}; 
        if contains(fileStr , 'x')
            tmp = strsplit(fileStr,{'x','ax'});
            expID = str2double(tmp{1});
            tmp2 = tmp{2};
            expID2 = str2double(tmp2(1:end-4));            
        else
            expID = str2double(fileStr(params.ID1s:end-params.ID1e));
            expID2 = str2double(fileStr(params.ID2s:end-params.ID2e));
        end
end


function [W] = ARMAsimple(ARMAtab , plotOn)
    M = 1;
    lag = 5;
    idList = [];
    for ii = 1:size(ARMAtab,1)-2
        if isnan(ARMAtab.rates(ii+1)) && ARMAtab.id2(ii) + 1 == ARMAtab.id2(ii+2)
            idList = [idList , ii+1];
        end
    end
    ARMAtab(idList,:) = [];
    samples = struct();
    samples.rates = (ARMAtab.rates - nanmean(ARMAtab.rates(:)))/nanstd(ARMAtab.rates(:));
    Mdl = varm(M , lag);
    [EstMdl,EstSE,logL,E] = estimate(Mdl , samples.rates);
    summary = summarize(EstMdl);
    W = reshape(cell2mat(EstMdl.AR),[M,M,lag]);
    if plotOn
        figure;
        bar(fliplr(reshape(W(1,1,:),[lag,1])'),'facecolor',rgb('red'))
        legend('amp -> amp','Location','northwest')
        xticks(1:lag)
        xlabel('time')
        xticklabels(-fliplr(1:lag))
    end
end
