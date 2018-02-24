tab = savTab;
params = getParams(DATASET); flags = getFlags(DATASET);

pathTo = params.TRACESpath;
fileList = getFileList(pathTo);
fixThresh = params.fixThresh;
% Siegel proposes to cutoff at 20%
LOWERBOUND = params.LOWERBOUND;
lAmps = [];
hAmps = [];
ARlist = []; pList = [];
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
                pList = [pList sTab.rates(kk)]; ARlist = [ARlist , nanmax(sMean)];
            end
        end
        pList = [pList nan]; ARlist = [ARlist , nan];
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
%%
ARMAplot(ARMAtab ,params, 1)
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

function [W] = ARMAplot(ARMAtab ,params, plotOn)
    M = 2; % number of features
    idList = [];
    for ii = 1:size(ARMAtab,1)-2
        if isnan(ARMAtab.rates(ii+1)) && ARMAtab.id2(ii) + 1 == ARMAtab.id2(ii+2)
            idList = [idList , ii+1];
        end
    end
    ARMAtab(idList,:) = [];


    samples = struct();
    %samples.rates = (ARMAtab.rates - nanmean(ARMAtab.rates(:)))/nanstd(ARMAtab.rates(:));
    samples.rates = (ARMAtab.rates - params.lEventUpper)/nanstd(ARMAtab.rates(:));
    samples.amps = (ARMAtab.amps - nanmean(ARMAtab.amps(:)))/nanstd(ARMAtab.amps(:));
    %samples.jitter = (ARMAtab.jitter - nanmean(ARMAtab.jitter(:)))/nanstd(ARMAtab.jitter(:));
    N = length(samples.rates);
    K = 7; % maximal number of lag
    ARs = 1:K;

    AICs = zeros(K,1); BICs = zeros(K,1); LogLs = zeros(K,1);
    for jj = 1:length(ARs)
        AR = ARs(jj);
        Mdl = varm(M , AR);
        EstMdl = estimate(Mdl , [samples.rates , samples.amps]);% , samples.jitter]);
        summary = summarize(EstMdl);
        AICs(jj) = summary.AIC; BICs(jj) = summary.BIC; LogLs(jj) = summary.LogLikelihood;
    end
    [mBI , mBID] = min(BICs);

    if plotOn
        figure;
        hold on
        plot([AICs , BICs])

        scatter(mBID ,  mBI , 'rx')
        xlabel('Autoregressive lag')
        ylabel('Information Criterion')
        legend('AIC','BIC')
        print(sprintf('fig_ARMA/%sAICBIC',''),'-dpng')
    end
    %mBID = 3*mBID;
    Mdl = varm(M , mBID);
    [EstMdl,EstSE,~,E] = estimate(Mdl , [samples.rates , samples.amps]);% , samples.jitter]);
    summary = summarize(EstMdl);
    %%
    W = reshape(cell2mat(EstMdl.AR),[M,M,mBID]);

    if plotOn
        figure; 
        subplot(1,2,1)
        hold on;
        err = reshape(table2array(summary.Table(3:end,2)),[M,M,mBID]);
        bar(fliplr(reshape(W(1,1,:),[mBID,1])'),'facecolor',rgb('orange'))
        bar(fliplr(reshape(W(1,2,:),[mBID,1])'),'facecolor',rgb('lightgrey'))
        errorbar(fliplr(reshape(W(1,2,:),[mBID,1])'),fliplr(reshape(err(1,2,:),[mBID,1])'),'x','color',rgb('black'))
        errorbar(fliplr(reshape(W(1,1,:),[mBID,1])'),fliplr(reshape(err(1,1,:),[mBID,1])'),'x','color',rgb('black'))
        ylim([-0.4 , 0.7])
        xticks(1:mBID)
        xticklabels(-fliplr(1:mBID))
        xlabel('time steps')
        ylabel('maximum likelihood coefficients')
        legend('amp -> rate','rate -> rate','Location','northwest')

        subplot(1,2,2)
        hold on;
        err = reshape(table2array(summary.Table(3:end,2)),[M,M,mBID]);
        bar(fliplr(reshape(W(2,2,:),[mBID,1])'),'facecolor',rgb('red'))
        bar(fliplr(reshape(W(2,1,:),[mBID,1])'),'facecolor',rgb('grey'))
        errorbar(fliplr(reshape(W(2,2,:),[mBID,1])'),fliplr(reshape(err(2,2,:),[mBID,1])'),'x','color',rgb('black'))
        errorbar(fliplr(reshape(W(2,1,:),[mBID,1])'),fliplr(reshape(err(2,1,:),[mBID,1])'),'x','color',rgb('black'))
        f1 = fit((1:mBID)' , fliplr(reshape(W(2,2,:),[mBID,1])')','exp1');
        plot(linspace(1,mBID+1),f1(linspace(1,mBID+1)),'--','color','black')
        f2 = fit((1:mBID)' , fliplr(reshape(W(2,1,:),[mBID,1])')','exp1');
        plot(linspace(1,mBID+1),f2(linspace(1,mBID+1)),'--','color','black')
        ylim([-0.4 , 0.7])
        legend('amp -> amp','rate -> amp','Location','northwest')
        xticks(1:mBID)
        xlabel('time')
        xticklabels(-fliplr(1:mBID))
    end

end
