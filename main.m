clear importCSV; 
close all
cd('/storage/gjor/Data/kirchnerj/src/DataAnalysis_final/')

addpath(genpath('tools'))
%%
% Wrapper for evaluating all regions
%
warning off
disp('Starting.')
DATASET = 'Fried';
params = getParams(DATASET); flags = getFlags(DATASET);
[savTab,savPcaTab,corTab,cpdTab] = combineFolder(params.TRACESpath,params , flags);
%%
params.COND = 0;
tab = savTab(logical((savTab.rates > params.lEventLower) + isnan(savTab.rates)) , :);
tab = tab(logical((tab.cond == params.COND) + (isnan(tab.cond))) , :); pcaTab = savPcaTab(savPcaTab.cond == params.COND , :);

%%
tab.amps = tab.amps + 1;


%%

uID = unique(tab.id(~isnan(tab.id)));

for ii = 1:length(uID)
    sTab = tab(tab.id == uID(ii),:);
    uID2 = unique(sTab.id2);
    for jj = 1:length(uID2)
        ssTab = sTab(sTab.id2 == uID2(jj),:);
        nF0 = nanmean(ssTab.amps(logical((ssTab.rates > 0.2).*(ssTab.rates < 0.8))));
        tab.amps(logical((tab.id == uID(ii)).*(tab.id2 == uID2(jj)))) = ssTab.amps/nF0;
    end
end

%%
scatterAmpJitter(tab,params,flags)% , data);
%%
boxParRate(tab,params)
%% 
dendoParAmp(tab,params,flags,'average','euclidean')
% 
%%
GMMcluster(tab,params,flags)
%%
% pause(1)
% numAnimalsBar(tab);
% pause(1)
% ageFrequencyBar(tab,params);
% pause(1)
% ampsBoxPlot(tab,params);
% pause(1)
% durationFWHMplot(tab,params);
% pause(1)
pcaPlot(pcaTab,params);
% pause(1)
% IEIs(tab , params);
% %% 
% if ~flags.skipBINOCULAR && ~flags.skipFORSKOLIN
%     pcaEnuFor(pcaTab,params);
% end
% correlationPlot(corTab,params);
% 
% correlationPooledPlot(corTab,params,flags);
%%
% params.lEventUpper = 0.8;
% LHrelationVaryDef(tab,params);
%% 
LHrelation(tab,params);
% 
% LHrelationPerID(tab,params);
% 
% LHrelationJitter(tab,params);

% tmp = LHgrid(savTab,params);

% tmp = LHgridSing(savTab,params);

% meanLplot(savTab,params, 0.2 , 0.8);

%%
W =  ARMAplot(tab ,params, flags.plotARMA);
% rec_W = shuffleARMA(tab , params);
%sigTestARMA(W,rec_W);
% 
% %%
%  ARMAsimple(tab, flags.plotARMA);
%%

% % cpdPlot(cpdTab);
f0plot(tab(~isnan(tab.id),:));
%%
cols = {rgb('red'),rgb('orange'),rgb('magenta'),rgb('pink')};
%  
% AgeFreqPlot3D(savTab,params,cols);
%  
 PCA3D(savPcaTab,params,cols);
% 
% Corr3D(corTab,params,cols);

%%

function [coeffmat] = LHgrid(savTab,params)
    N = 15;
    coeffmat = zeros(N,N);
    dl = linspace(0.2,0.55,N); du = linspace(0.6 , 0.9,N);
    for ll = 1:length(dl)
        for uu = 1:length(du)
            tab = savTab;
            uID = unique(tab.id(~isnan(tab.id)));
            for ii = 1:length(uID)
                sTab = tab(tab.id == uID(ii),:);
                nF0 = nanmean(sTab.amps(logical((sTab.rates > dl(ll)).*(sTab.rates < du(uu)))));
                tab.amps(tab.id == uID(ii)) = sTab.amps/nF0;

                uID2 = unique(sTab.id2);
                for jj = 1:length(uID2)
                    ssTab = sTab(sTab.id2 == uID2(jj),:);
                    nF0 = nanmean(ssTab.amps(logical((ssTab.rates > dl(ll)).*(ssTab.rates < du(uu)))));
                    if isnan(nF0)
                        tab(logical((tab.id == uID(ii)).*(tab.id2 == uID2(jj))),:) = [];
                    else
                        tab.amps(logical((tab.id == uID(ii)).*(tab.id2 == uID2(jj)))) = ssTab.amps/nF0;                        
                    end
                end
            end
            tab.amps = tab.amps + 1;
            coeffmat(ll,uu) = LHrelation(tab,params,0);
        end
    end
    figure();
    imagesc(coeffmat');
    xticks(1:3:N);yticks(1:3:N);
    xticklabels(round(du(1:3:N),2)); yticklabels(round(dl(1:3:N),2))
    colorbar;
end


function [coeffmat] = LHgridSing(savTab,params)
    N = 15;
    coeffmat = zeros(N,N);
    dl = linspace(0.2,0.55,N); du = linspace(0.6 , 0.9,N);
    for ll = 1:length(dl)
        for uu = 1:length(du)
            tab = savTab;
            uID = unique(tab.id(~isnan(tab.id)));
            for ii = 1:length(uID)
                sTab = tab(tab.id == uID(ii),:);
                nF0 = nanmean(sTab.amps(logical((sTab.rates > dl(ll)).*(sTab.rates < du(uu)))));
                tab.amps(tab.id == uID(ii)) = sTab.amps/nF0;

                uID2 = unique(sTab.id2);
                for jj = 1:length(uID2)
                    ssTab = sTab(sTab.id2 == uID2(jj),:);
                    nF0 = nanmean(ssTab.amps(logical((ssTab.rates > dl(ll)).*(ssTab.rates < du(uu)))));
                    if isnan(nF0)
                        tab(logical((tab.id == uID(ii)).*(tab.id2 == uID2(jj))),:) = [];
                    else
                        tab.amps(logical((tab.id == uID(ii)).*(tab.id2 == uID2(jj)))) = ssTab.amps/nF0;                        
                    end
                end
            end
            tab.amps = tab.amps + 1;
            coeffmat(ll,uu) = LHrelationSing(tab,params,0);
        end
    end
    figure();
    imagesc(coeffmat');
    xticks(1:3:N);yticks(1:3:N);
    xticklabels(round(du(1:3:N),2)); yticklabels(round(dl(1:3:N),2))
    colorbar;
end

function ccoeffs = LHrelationVaryDef(tab,params)
    N = 100;
    ccoeffs = zeros(N,1);
    dx = linspace(0.3,0.99,N);
    for ii = 1:N
        params.lEventUpper = dx(ii);
        ccoeffs(ii) = LHrelation(tab,params,0);
    end
    figure; plot(dx,ccoeffs,'LineWidth',4)
    xlabel('H Event Threshold')
    ylabel('Correlation Coefficient')
    title('Correlation Coefficient for different H-event Thresholds')
    %ylim([0,1])
end

function [] = GMMcluster(tab,params,flags)
    X = [tab.amps + 1 , tab.jitter , tab.rates];
    options = statset('Display','final');
    obj = fitgmdist(X,1,'Replicates',20,'Options',options);
    ppb = posterior(obj , [tab.amps + 1 , tab.jitter , tab.rates]);
    [~,idx] = max(ppb');
    figure;
    subplot(2,1,1)
    colormap(jet); scatter3(tab.amps + 1 , tab.jitter , tab.rates ,20,rgb('darkblue'), 'filled')
    xlabel('Amplitude (F/F0)')
    ylabel(sprintf('Jitter \n'))
    zlabel('Participation rate')
    title(sprintf('AIC: %4.0f , BIC: %4.0f',obj.AIC , obj.BIC))

    obj = fitgmdist(X,3,'Replicates',20,'Options',options);
    ppb = posterior(obj , [tab.amps + 1 , tab.jitter , tab.rates]);
    [~,idx] = max(ppb');
    subplot(2,1,2)
    colormap(jet); scatter3(tab.amps + 1 , tab.jitter , tab.rates ,20, idx, 'filled')
    xl = xlim(); yl = ylim(); zl = zlim();
    hold on
    f = fsurf(0.6);
    alpha(f,0.25)
    f = fsurf(0.8);
    alpha(f,0.25)
    title(sprintf('AIC: %4.0f , BIC: %4.0f',obj.AIC , obj.BIC))
    xlim(xl);ylim(yl);zlim(zl);

end

function [] = Corr3D(corTab,params,cols )
    addpath('tools')
    uCond = unique(corTab.cond(~isnan(corTab.cond)));
    uAge = unique(corTab.age(~isnan(corTab.age)));
    uDist = unique(corTab.range(~isnan(corTab.range)));
    corMat = zeros(length(uCond),length(uAge),length(uDist));
    for ii = 1:length(uCond)
        cCond = uCond(ii);
        for jj = 1:length(uAge)
            cAge = uAge(jj);
            for kk = 1:length(uDist)
                cDist = uDist(kk);
                corMat(ii,jj,kk) = nanmean(corTab.cors(logical((corTab.age == cAge) .* (corTab.cond == cCond) .* (corTab.range == cDist))));
            end
        end
    end
    corMat(corMat == 0) = NaN;
    rcorMat = nanmean(corMat,3);
    figure;subplot(1,2,1);
    h = bar3(rcorMat'); 
    for ii = 1:length(uCond)
        set(h(ii),'FaceColor',cols{ii});
    end
    xlim([0.5 , length(uCond) + 0.5]); ylim([0.5 , length(uAge) + 0.5]);
    xticklabels(params.labels); yticklabels(cellstr(num2str(uAge)));
    xlabel('Area'); ylabel('Postnatal day'); zlabel('Mean Correlation Coefficient')
    subplot(1,2,2);hold on;
    axs = [];
    for ii = 1:length(uCond)
        for jj = 1:length(uDist)
            axs(ii) = plot3(1:length(uAge),repmat(jj,length(uAge),1),corMat(ii,:,jj),'-o','Color',cols{ii},'LineWidth',3);
        end
    end
    xticks(1:length(uAge));xticklabels(cellstr(num2str(uAge)));xlabel('Postnatal day')
    yticks(1:length(uDist));
    yticklabels({sprintf('Corr %d-%d \\mu m', params.lNeigh , params.mNeigh),...
        sprintf('Corr %d-%d \\mu m', params.mNeigh , params.umNeigh),...
        sprintf('Corr %d-%d \\mu m', params.umNeigh , params.uNeigh),...
        sprintf('Corr > %d \\mu m', params.uNeigh)})
    legend(axs,params.labels,'Location','Northeast');ylabel('Distance');zlabel('Mean Correlation Coefficient')
    view(3)
end

function [] = AgeFreqPlot3D(tab,params,cols )
    tab(tab.rates < params.lEventLower , : ) = [];
    addpath('tools')
    uCond = unique(tab.cond(~isnan(tab.cond)));
    uAge = unique(tab.age(~isnan(tab.age)));

    freqAgeCond = zeros(length(uCond),length(uAge));
    for ii = 1:length(uCond)
        cCond = uCond(ii);
        [cAges , relMult] = ageFrequencyBar(tab(logical((tab.cond == cCond) + isnan(tab.cond)),:),params,0);
        for jj = 1:length(uAge)
            if ~isempty(cAges) && uAge(jj) == cAges(1)
                freqAgeCond(ii , jj) = relMult(1);
                cAges(1) = []; relMult(1) = [];
            end
        end
    end
    freqAgeCond(freqAgeCond == 0) = NaN;
    figure;subplot(1,2,1);
    h = bar3(freqAgeCond'); 
    for ii = 1:length(uCond)
        set(h(ii),'FaceColor',cols{ii});
    end
    xlim([0.5 , length(uCond) + 0.5]); ylim([0.5 , length(uAge) + 0.5]);
    xticklabels(params.labels); yticklabels(cellstr(num2str(uAge)));
    xlabel('Area'); ylabel('Postnatal day'); zlabel('Frequency (1/min)')
    subplot(1,2,2);hold on;
    for ii = 1:length(uCond)
        plot(1:length(uAge),freqAgeCond(ii,:),'-o','Color',cols{ii},'LineWidth',3)
    end
    xticks(1:length(uAge));xticklabels(cellstr(num2str(uAge)));xlabel('Postnatal day')
    legend(params.labels,'Location','Northwest');ylabel('Frequency (1/min)')
end

function [] = PCA3D(pcaTab,params,cols)
    uCond = unique(pcaTab.cond(~isnan(pcaTab.cond)));
    uAge = unique(pcaTab.age(~isnan(pcaTab.age)));

    PCAs = zeros(length(uCond),length(uAge));
    for ii = 1:length(uCond)
        cCond = uCond(ii);
        for jj = 1:length(uAge)
            cAge = uAge(jj);
            PCAs(ii,jj) = nanmean(pcaTab.pcaPerc(logical((pcaTab.age == cAge) .* (pcaTab.cond == cCond))));
        end
    end
    PCAs(PCAs == 0) = NaN;
    figure;subplot(1,2,1);
    h = bar3(PCAs'); xlim([0.5 , length(uCond) + 0.5]); ylim([0.5 , length(uAge) + 0.5]);
    for ii = 1:length(uCond)
        set(h(ii),'FaceColor',cols{ii});
    end
    xticklabels(params.labels); yticklabels(cellstr(num2str(uAge)));
    xlabel('Area'); ylabel('Postnatal day'); zlabel('Average dimensionality')
    subplot(1,2,2);hold on;
    for ii = 1:length(uCond)
        plot(1:length(uAge),PCAs(ii,:),'-o','Color',cols{ii},'LineWidth',3)
    end
    xticks(1:length(uAge));xticklabels(cellstr(num2str(uAge)));xlabel('Postnatal day')
    legend(params.labels,'Location','Northwest');ylabel('Average dimensionality')
end

function [] = dendoParAmp(tab,params,flags,method,method2)
    figure;
    labs = floor(tab.rates*9.99);
    mAmps = zeros(10,1); mJitter = zeros(10,1); mPar = zeros(10,1);
    for ii = 0:9
        mAmps(ii+1) = nanmean(tab.amps(labs == ii));
        mJitter(ii+1) = nanmean(tab.jitter(labs == ii));
        mPar(ii+1) = nanmean(tab.rates(labs == ii));
    end
    L = linkage([mAmps mJitter],method,method2); dendrogram(L,'Reorder',1:10);
    xticklabels({'0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%'});
    title(sprintf('M1: %s , M2: %s',method , method2))
end

function [] = meanLplot(tab,params, lNorm , hNorm)
    N = 15;
    tab.amps = tab.amps + 1;
    uID = unique(tab.id);
    uID = uID(~isnan(uID));
    %cm = jet(5);
    figure;hold on; 
    for ii = 1:N%length(uID)
        sTab = tab(tab.id == uID(ii),:);
        uID2 = unique(sTab.id2);
        uID2 = uID2(~isnan(uID2));
        recAmps = [];
        recAmps2 = [];

        Ts = [0];
        for jj = 1:length(uID2)
            ssTab = sTab(sTab.id2 == uID2(jj) , :);
            recAmps = [recAmps nanmean(ssTab.amps(logical((ssTab.rates > lNorm).*(ssTab.rates < hNorm))))];
            recAmps2 = [recAmps2 nanmean(ssTab.amps(logical((ssTab.rates > hNorm))))];

            Ts = [Ts , Ts(end)+ ssTab.T(1)];
        end
        Ts = Ts(2:end);
        plot(Ts , recAmps,'-x');%,'Color',cm(ii,:))
        %plot(Ts + 50 , recAmps2,'-x','Color',cm(ii,:))

        %plot([0 , 35],[ii ii],'Color',rgb('gray'))
    end
end

function [] = f0plot(tab)
    uID = unique(tab.id);
    figure;hold on; 
    for ii = 1:length(uID)
        cF0 = unique(tab.F0quant(tab.id == uID(ii)));
        uID2 = unique(tab.id2(tab.id == uID(ii)));

        plot(uID2 , (cF0-min(cF0))/max(cF0-min(cF0))+ii,'-x','Color','black')
        plot([0 , 18],[ii ii],'Color',rgb('gray'))
    end
end

function [] = sigTestARMA(W,rec_W)
    figure; c = 1;
    for ii = 1:size(W,3)
        for jj = 1:size(W,2)
           for kk = 1:size(W,1)
               subplot(4,4,c)
               samp = rec_W(kk,jj,ii,:);
               samp = samp(:);
               hold on;title(sprintf('lag: %d , W(%d,%d)',ii,kk,jj))
               hist(samp); plot([W(kk,jj,ii) , W(kk,jj,ii)] , [0,100],'r','LineWidth',5)
               c = c + 1;
           end
        end
    end

end

function [rec_W] = shuffleARMA(tab,params)
    M = params.shuffleNum;
    removeIDs = [];

    for ii = 2:height(tab)-1
        contCond = isnan(tab.amps(ii)) && tab.id(ii+1) == tab.id(ii-1) ...
            && tab.id2(ii+1) == tab.id2(ii-1) + 1;
        if contCond
           removeIDs = [removeIDs , ii]; 
        end
    end
    localTab = tab; localTab(removeIDs,:) = [];
    uNaN = [];
    for jj = 2:height(localTab)-1
        contCond = isnan(localTab.amps(jj)) && localTab.id(jj+1) ~= localTab.id(jj-1);
        if contCond
           uNaN = [uNaN , jj]; 
        end
    end
    rec_W = zeros(2,2,4,M);
    for ii = 1:M
        ii
        shuffleTab = localTab;
        for jj =2:length(uNaN)
            pID = uNaN(jj-1); nID = uNaN(jj);
            cRow = shuffleTab(pID+1:nID-1 , :);
            shuffleTab(pID+1:nID-1 , :) = cRow(randperm(height(cRow)),:);
        end
        [W] = ARMAplot(shuffleTab ,params, 0);
        rec_W(:,:,1:min(size(W,3),4),ii) = W(:,:,1:min(size(W,3),4));
    end
end

function [] = boxParRate(tab,params)
    figure();subplot(1,2,1);hold on
    IDXs = zeros(length(tab.amps),1);
    boxplot(tab.amps + 1 , floor(tab.rates*(4.999)));
    xticklabels({'0-20%','20-40%','40-60%','60-80%','80-100%'});
    xlabel('Participation rate')
    ylabel('Amplitude (F/F0)')
%     for ii = 1:4
%         [h,p,ci,stats] = ttest2(tab.amps(floor(tab.rates*(5.999)) == ii),tab.amps(floor(tab.rates*(5.999)) == ii+1));
%         if h == 1        
%            sigstar({[ii,ii+1]},p*4)
%         end
%     end
    subplot(1,2,2)
    boxplot(tab.jitter, floor(tab.rates*(4.999)));
    xticklabels({'0-20%','20-40%','40-60%','60-80%','80-100%'});
    xlabel('Participation rate')
    ylabel('Jitter (s)')
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
    samples.jitter = (ARMAtab.jitter - nanmean(ARMAtab.jitter(:)))/nanstd(ARMAtab.jitter(:));
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

function LHrelationJitter(LHtab , params)
    removeIDs = [];

    for ii = 2:height(LHtab)-1
        contCond = isnan(LHtab.amps(ii)) && LHtab.id(ii+1) == LHtab.id(ii-1) ...
            && LHtab.id2(ii+1) == LHtab.id2(ii-1) + 1;
        if contCond
           removeIDs = [removeIDs , ii]; 
        end
    end
    rLHtab = LHtab; rLHtab(removeIDs,:) = [];
    %%

    hIDs = find(rLHtab.rates > params.lEventUpper);
    hJit = rLHtab.jitter(hIDs) + 1; hJit = hJit(2:end);
    lJit = zeros(length(hIDs)-1 , 1);
    IDs = zeros(length(hIDs)-1 , 1);
    for ii = 2:length(hIDs)
        prevID = hIDs(ii-1) + 1;
        nextID = hIDs(ii) - 1;
        IDs(ii-1) = rLHtab.id(nextID);
        slice = rLHtab.jitter(prevID:nextID);
        if isempty(slice)
            lJit(ii - 1) = NaN;
        else
            bID = find(isnan(slice));
            if isempty(bID)
                lJit(ii - 1) = mean(slice);
            else
                bID = bID(end);
                lJit(ii - 1) = mean(slice(bID+1:end));
            end
        end
    end
    lJit = lJit + 1;
    %%
    rlJitter = lJit(logical((~isnan(lJit)).*(~isnan(hJit))));
    rhJitter = hJit(logical((~isnan(lJit)).*(~isnan(hJit))));
    rIDs = IDs(logical((~isnan(lJit)).*(~isnan(hJit))));
    figure()
    scatter(rlJitter,rhJitter ,45,rIDs,'filled')
    lsline
    [R,P,RL,RU] = corrcoef(rlJitter,rhJitter)
    xlabel('Mean L-event Jitter (s)')
    ylabel('H-event Jitter (s)')

end

function [] =  IEIs(tab , params)
    %%
    removeIDs = [];

    for ii = 2:height(tab)-1
        contCond = isnan(tab.amps(ii)) && tab.id(ii+1) == tab.id(ii-1) ...
            && tab.id2(ii+1) == tab.id2(ii-1) + 1;
        if contCond
           removeIDs = [removeIDs , ii]; 
        end
    end
    LHtab = tab; LHtab(removeIDs,:) = [];
    figure();hold on;
    dxx = linspace(0,300,31);
    HEvent = double(LHtab.rates > params.lEventUpper);

    HIEI = diff(LHtab.realStartTimes(find(HEvent)));
    HIEI = HIEI(HIEI > 0);
    [NHIEI,~] = histcounts(HIEI,dxx);
    LIEI = diff(LHtab.realStartTimes(find(1 - HEvent)));
    LIEI = LIEI(~isnan(LIEI));
    [NLIEI,~] = histcounts(LIEI,dxx);
    dx = linspace(0,300,60);
    Lbins = zeros(60,1);    Hbins = zeros(60,1);
    Lbins(1:2:end) = NLIEI;
    Hbins(2:2:end) = NHIEI;

    histogram('BinEdges',dx,'BinCounts',Lbins(1:end-1),'Normalization','Probability','FaceColor',rgb('darkgray'))
    histogram('BinEdges',dx,'BinCounts',Hbins(1:end-1),'Normalization','Probability','FaceColor',rgb('red'))
    xlim([0,300])
    legend({'Participation 20-80%','Participation 80-100%'})
    legend('boxoff')
end

function [ccoeff] = LHrelation(LHtab,params,plotFlag)
    if nargin < 3
        plotFlag = 1;
    end
    removeIDs = [];

    for ii = 2:height(LHtab)-1
        contCond = isnan(LHtab.amps(ii)) && LHtab.id(ii+1) == LHtab.id(ii-1) ...
            && LHtab.id2(ii+1) == LHtab.id2(ii-1) + 1;
        if contCond
           removeIDs = [removeIDs , ii]; 
        end
    end
    rLHtab = LHtab; rLHtab(removeIDs,:) = [];
    %%

    hIDs = find(rLHtab.rates > params.lEventUpper);
    hAmp = rLHtab.amps(hIDs); hAmp = hAmp(2:end);
    lAmp = zeros(length(hIDs)-1 , 1);
    IDs = zeros(length(hIDs)-1 , 1);
    AGEs = zeros(length(hIDs)-1 , 1); 
    for ii = 2:length(hIDs)
        prevID = hIDs(ii-1) + 1;
        nextID = hIDs(ii) - 1;
        IDs(ii-1) = rLHtab.id(nextID);
        AGEs(ii-1) = rLHtab.age(nextID);
        slice = rLHtab.amps(prevID:nextID);
        if isempty(slice)
            lAmp(ii - 1) = NaN;
        else
            bID = find(isnan(slice));
            if isempty(bID)
                lAmp(ii - 1) = mean(slice);
            else
                bID = bID(end);
                lAmp(ii - 1) = mean(slice(bID+1:end));
            end
        end
    end
    %%
    rlAmp = lAmp(~isnan(lAmp));
    rhAmp = hAmp(~isnan(lAmp));
    rIDs = IDs(~isnan(lAmp));
    rAGES = AGEs(~isnan(lAmp));
    
    [R,P,RL,RU] = corrcoef(rlAmp,rhAmp);
    ccoeff = R(1,2);
    if plotFlag
        figure()
        scatter(rlAmp,rhAmp ,45,'MarkerFaceColor',rgb('gray'),'MarkerEdgeColor',rgb('black'))
        lsline
        %colorbar
        title(sprintf('H-event: %2.0f , correlation: %1.2f',params.lEventUpper, ccoeff))
        xlabel('Mean L-event amplitude ^{F}/_{F0}')
        ylabel('H-event amplitude ^{F}/_{F0}')
    end
end

function [ccoeff] = LHrelationSing(LHtab,params,plotFlag)
    if nargin < 3
        plotFlag = 1;
    end
    removeIDs = [];

    for ii = 2:height(LHtab)-1
        contCond = isnan(LHtab.amps(ii)) && LHtab.id(ii+1) == LHtab.id(ii-1) ...
            && LHtab.id2(ii+1) == LHtab.id2(ii-1) + 1;
        if contCond
           removeIDs = [removeIDs , ii]; 
        end
    end
    rLHtab = LHtab; rLHtab(removeIDs,:) = [];
    %%

    hIDs = find(rLHtab.rates > params.lEventUpper);
    hAmp = rLHtab.amps(hIDs); hAmp = hAmp(2:end);
    lAmp = zeros(length(hIDs)-1 , 1);
    IDs = zeros(length(hIDs)-1 , 1);
    AGEs = zeros(length(hIDs)-1 , 1); 
    for ii = 2:length(hIDs)
        prevID = hIDs(ii-1) + 1;
        nextID = hIDs(ii) - 1;
        IDs(ii-1) = rLHtab.id(nextID);
        AGEs(ii-1) = rLHtab.age(nextID);
        slice = rLHtab.amps(prevID:nextID);
        if isempty(slice)
            lAmp(ii - 1) = NaN;
        else
            bID = find(isnan(slice));
            if isempty(bID)
                lAmp(ii - 1) = mean(slice(end));
            else
                bID = bID(end);
                lAmp(ii - 1) = mean(slice(end));
            end
        end
    end
    %%
    rlAmp = lAmp(~isnan(lAmp));
    rhAmp = hAmp(~isnan(lAmp));
    rIDs = IDs(~isnan(lAmp));
    rAGES = AGEs(~isnan(lAmp));
    
    [R,P,RL,RU] = corrcoef(rlAmp,rhAmp);
    ccoeff = R(1,2);
    if plotFlag
        figure()
        scatter(rlAmp,rhAmp ,45,'MarkerFaceColor',rgb('gray'),'MarkerEdgeColor',rgb('black'))
        lsline
        %colorbar
        title(sprintf('H-event: %2.0f , correlation: %1.2f',params.lEventUpper, ccoeff))
        xlabel('Mean L-event amplitude ^{F}/_{F0}')
        ylabel('H-event amplitude ^{F}/_{F0}')
    end
end

function [] = LHrelationPerID(LHtab,params)
    removeIDs = [];

    for ii = 2:height(LHtab)-1
        contCond = isnan(LHtab.amps(ii)) && LHtab.id(ii+1) == LHtab.id(ii-1) ...
            && LHtab.id2(ii+1) == LHtab.id2(ii-1) + 1;
        if contCond
           removeIDs = [removeIDs , ii]; 
        end
    end
    rLHtab = LHtab; rLHtab(removeIDs,:) = [];
    %%

    hIDs = find(rLHtab.rates > params.lEventUpper);
    hAmp = rLHtab.amps(hIDs); hAmp = hAmp(2:end);
    lAmp = zeros(length(hIDs)-1 , 1);
    IDs = zeros(length(hIDs)-1 , 1);
    AGEs = zeros(length(hIDs)-1 , 1);
    F0s = zeros(length(hIDs)-1 , 1);
    for ii = 2:length(hIDs)
        prevID = hIDs(ii-1) + 1;
        nextID = hIDs(ii) - 1;
        IDs(ii-1) = rLHtab.id(nextID);
        AGEs(ii-1) = rLHtab.age(nextID);
        F0s(ii-1) = rLHtab.F0(nextID);

        slice = rLHtab.amps(prevID:nextID);
        if isempty(slice)
            lAmp(ii - 1) = NaN;
        else
            bID = find(isnan(slice));
            if isempty(bID)
                lAmp(ii - 1) = mean(slice);
            else
                bID = bID(end);
                lAmp(ii - 1) =mean(slice(bID+1:end));
            end
        end
    end
    %%
    rlAmp = lAmp(~isnan(lAmp));
    rhAmp = hAmp(~isnan(lAmp));
    rIDs = IDs(~isnan(lAmp));
    rAGEs = AGEs(~isnan(lAmp));
    rF0s = F0s(~isnan(lAmp));

    uIDs = unique(rIDs);
    figure(); hold on;
    corrs = [];
    markers = {'o'};%,'*','.','x','s','d','^','v','>','<','p','h'};
    for ii = 1:length(uIDs)
        scatter(rlAmp(rIDs == uIDs(ii)),rhAmp(rIDs == uIDs(ii)) ,45 , markers{mod(ii,numel(markers))+1} ,'filled')
        %lsline
        corrs = [corrs , corr(rlAmp(rIDs == uIDs(ii)),rhAmp(rIDs == uIDs(ii)) )];
    end
%    nanmean(corrs)
%    figure; histogram(corrs)
%     subplot(4,ceil(length(uIDs)/4),ii+1); hold on;
%     scatter(rlAmp,rhAmp ,45,rIDs,'filled')
%     lsline
    [R,P,RL,RU] = corrcoef(rlAmp,rhAmp)
    [b,bint] = regress(rhAmp, [ones(length(rlAmp),1), rlAmp , rAGEs , rF0s])
    xlabel('Mean L-event amplitude ^{F}/_{F0}')
    ylabel('H-event amplitude ^{F}/_{F0}')
end

function [] = pcaEnuFor(pcaTab , params)
    figure;
    subplot(1,2,1)
    hold on;
    enuTab = pcaTab(pcaTab.cond == 1,:);
    enuID = unique(enuTab.id);
    bEnu = [];
    aEnu = [];
    for ii = 1:length(enuID)
        cID = enuID(ii);
        bEnu = [bEnu , nanmean(pcaTab.pcaPerc(logical((pcaTab.id == cID) .* (pcaTab.cond == 0))))];
        aEnu = [aEnu , nanmean(pcaTab.pcaPerc(logical((pcaTab.id == cID) .* (pcaTab.cond == 1))))];
        plot([1,2],[bEnu(end),aEnu(end)],'color',rgb('black'))
        %fullTab = [fullTab ; pcaTab(pcaTab.id == cID , :)];
    end
    scatter(ones(length(bEnu),1),bEnu,'filled','MarkerEdgeColor',rgb('black'))
    scatter(2*ones(length(bEnu),1),aEnu,'filled','MarkerEdgeColor',rgb('black'))
    xlim([0,3]);ylim([0.15,0.4]);ylabel('% of full dimensionality')
    xticks([1,2]);xticklabels({'before EN','after EN'});xtickangle(45)
    [~,p] = ttest(bEnu,aEnu)
    sigstar({[1,2]},p);

    subplot(1,2,2)
    hold on;
    forTab = pcaTab(pcaTab.cond == 2,:);
    forID = unique(forTab.id);
    bFor = [];
    aFor = [];
    for ii = 1:length(forID)
        cID = forID(ii);
        bFor = [bFor , nanmean(pcaTab.pcaPerc(logical((pcaTab.id == cID) .* (pcaTab.cond == 0))))];
        aFor = [aFor , nanmean(pcaTab.pcaPerc(logical((pcaTab.id == cID) .* (pcaTab.cond == 2))))];
        plot([1,2],[bFor(end),aFor(end)],'color',rgb('black'))
        %fullTab = [fullTab ; pcaTab(pcaTab.id == cID , :)];
    end
    scatter(ones(length(bFor),1),bFor,'filled','MarkerEdgeColor',rgb('black'))
    scatter(2*ones(length(bFor),1),aFor,'filled','MarkerEdgeColor',rgb('black'))
    xlim([0,3]);ylim([0.15,0.4])
    xticks([1,2]);xticklabels({'before FO','after FO'});xtickangle(45);
    [~,p] = ttest(bFor,aFor)
    sigstar({[1,2]},p);
end

function [] = scatterAmpJitter(tab,params,flags , data)
    figure(); 
    subplot(1,2,1);
    colormap(jet);scatter(tab.rates,tab.amps,params.ScatterPointSize,tab.durationsFWHM,'filled')%,'MarkerFaceColor',rgb('gray'),'MarkerEdgeColor',rgb('black')); 
    colorbar
    ylabel('Amplitude ^{F}/_{F0}'); xlim([0,1]); %ylim([0.8,1.6])
    xlabel('Participation rate (%)');% colorbar
    
    if nargin > 3
        hold on 
        scatter(data(:,7)/100,data(:,6),20,'red');%,'filled')
    end
        
    if flags.addRegLine
        hold on;
        X = [ones(length(tab.rates(tab.rates < params.lEventUpper)),1) tab.rates(tab.rates < params.lEventUpper)];
        b = X\(tab.amps(tab.rates < params.lEventUpper)+1);
        dp = linspace(params.lEventLower,params.lEventUpper);
        plot(dp , dp*b(2) + b(1),'r','LineWidth',2)
        X = [ones(length(tab.rates(tab.rates > params.lEventUpper)),1) tab.rates(tab.rates > params.lEventUpper)];
        b = X\(tab.amps(tab.rates > params.lEventUpper)+1);
        dp = linspace(params.lEventUpper,1);
        plot(dp , dp*b(2) + b(1),'b','LineWidth',2)
    end
    
    if flags.addRegLine
        hold on;
        X = [ones(length(partvsamp.Participation(partvsamp.Participation < params.lEventUpper)),1) partvsamp.Participation(partvsamp.Participation < params.lEventUpper)];
        b = X\(partvsamp.FF0(partvsamp.Participation < params.lEventUpper));
        dp = linspace(params.lEventLower,params.lEventUpper);
        plot(dp , dp*b(2) + b(1),'r','LineWidth',2)
        X = [ones(length(partvsamp.Participation(partvsamp.Participation > params.lEventUpper)),1) partvsamp.Participation(partvsamp.Participation > params.lEventUpper)];
        b = X\(partvsamp.FF0(partvsamp.Participation > params.lEventUpper));
        dp = linspace(params.lEventUpper,1);
        plot(dp , dp*b(2) + b(1),'b','LineWidth',2)
    end
    subplot(1,2,2);
    colormap(jet);scatter(tab.rates,tab.jitter,params.ScatterPointSize,'MarkerFaceColor',rgb('gray'),'MarkerEdgeColor',rgb('black')); 
    xlabel('Participation Rate'); ylabel('Jitter'); xlim([0,1]);
    if nargin > 3
        hold on 
        scatter(data(:,7)/100,data(:,5),20,'red');%,'filled')
    end
%     print(sprintf('fig_Fried/%sScatter','V1'),'-dpng')

end

function [] = numAnimalsBar(tab)
    figure();
    uniqueAges = unique(tab.age(~isnan(tab.age)));
    numA = zeros(length(uniqueAges),1);
    for ii = 1:length(uniqueAges)
        uAge = uniqueAges(ii);
        uID = unique(tab.id(tab.age == uAge));
        numA(ii) =  length(uID);
    end
    bar(numA)
    xticklabels(uniqueAges)
    xlabel('Age of animal in days')
    ylabel('Count')
    title('Number of animals recorded')
end

function [uAges , relMult] =  ageFrequencyBar(tab,params,plotOn)
    if nargin < 3
        plotOn = 1;
    end
    if plotOn
        figure();hold on;
    end
    uAges = unique(tab.age(~isnan(tab.age)));
    
    mult = histc(tab.age , uAges);
    hMult = histc(tab.age(tab.rates > params.lEventUpper) , uAges);
    lMult = histc(tab.age(tab.rates < params.lEventUpper) , uAges);
    relMult = zeros(length(mult),1);
    relHMult = zeros(length(hMult),1);
    relLMult = zeros(length(lMult),1);

    for ii = 1:length(uAges)
        totalT = 0;
        cAge = uAges(ii);
        aTab = tab(tab.age == cAge , :);
        uID = unique(aTab.id);
        for jj = 1:length(uID)
            cID = uID(jj);
            uID2 = unique(aTab.id2(aTab.id == cID));
            for kk = 1:length(uID2)
                ccID = uID2(kk);
                totalT = totalT + unique(aTab.T(logical((aTab.id == cID) .* (aTab.id2 == ccID))))/60;
            end    
        end
        
        relMult(ii) = mult(ii)/totalT; 
        relHMult(ii) = hMult(ii)/totalT; 
        relLMult(ii) = lMult(ii)/totalT; 
    end
    if plotOn
        bar(uAges , relMult)
        bar(uAges , relLMult)
        bar(uAges , relHMult)
        xlabel('Postnatal days')
        ylabel('frequency')
        legend({'pooled','L-events','H-event'},'Location','Northwest')
    end
end

function [] = ampsBoxPlot(tab,params)   
    figure();
    subplot(1,3,1)
    boxplot(tab.amps(~isnan(tab.amps)) + 1,tab.age(~isnan(tab.amps)))
    ylabel('Amplitude')
    title('# all events > 20%')

    subplot(1,3,2)

    boxplot(tab.amps(tab.rates <= params.lEventUpper) + 1,tab.age(tab.rates <= params.lEventUpper))
    xlabel('Postnatal days')
    title('# L events \leq 80%')

    subplot(1,3,3)

    boxplot(tab.amps(tab.rates > params.lEventUpper) + 1,tab.age(tab.rates > params.lEventUpper))
    title('# H events > 80%')

    print(sprintf('fig_Fried/%sAmps','V1'),'-dpng')
end

function [] = durationFWHMplot(tab,params)
    figure();
    title('Durations of events, FWHM')
    subplot(1,3,1)
    boxplot(tab.durationsFWHM,tab.age)
    ylabel('Duration')
    title('# all events > 20%')
    ylim([0,14])

    subplot(1,3,2)

    boxplot(tab.durationsFWHM(tab.rates <= params.lEventUpper),tab.age(tab.rates <= params.lEventUpper))
    xlabel('Postnatal days')
    title('# L events \leq 80%')
    ylim([0,14])

    subplot(1,3,3)

    boxplot(tab.durationsFWHM(tab.rates > params.lEventUpper),tab.age(tab.rates > params.lEventUpper))
    title('# H events > 80%')
    ylim([0,14])
    print(sprintf('fig_Fried/%sDurations','V1'),'-dpng')
end

function [] = pcaPlot(pcaTab,params)
    figure()
    boxplot(pcaTab.pcaPercAlt,pcaTab.age)
    title(sprintf('Data dimensionality as the percentage\n of PCs required to explain 80% of variance'))
    xlabel('Postnatal days')
    ylabel('% of full dimensionality')
    print(sprintf('fig_Fried/%sPCA','V1'),'-dpng')
end

function [] = correlationPlot(corTab,params)   
    figure(); 
    hold on
    title('Short range pairwise correlations')
    boxplot(corTab.cors(corTab.range == 1),corTab.age(corTab.range == 1))
    xlabel('Postnatal days')
    ylabel(sprintf('Corr %d-%d \\mu m', params.lNeigh , params.mNeigh))
    print(sprintf('fig_Fried/%sCorrs1','V1'),'-dpng')
    figure(); 
    hold on
    title('Mid range pairwise correlations')
    boxplot(corTab.cors(corTab.range == 2),corTab.age(corTab.range == 2))
    xlabel('Postnatal days')
    ylabel(sprintf('Corr %d-%d \\mu m', params.mNeigh , params.umNeigh))
    print(sprintf('fig_Fried/%sCorrs2','V1'),'-dpng')
    figure(); 
    hold on
    title('Long range pairwise correlations')
    boxplot(corTab.cors(corTab.range == 3),corTab.age(corTab.range == 3))
    xlabel('Postnatal days')
    ylabel(sprintf('Corr %d-%d \\mu m', params.umNeigh , params.uNeigh))
    print(sprintf('fig_Fried/%sCorrs3','V1'),'-dpng')
    figure(); 
    hold on
    title('Very Long range pairwise correlations')
    boxplot(corTab.cors(corTab.range == 4),corTab.age(corTab.range == 4))
    xlabel('Postnatal days')
    ylabel(sprintf('Corr > %d \\mu m', params.uNeigh))
    print(sprintf('fig_Fried/%sCorrs4','V1'),'-dpng')
end

function [] = correlationPooledPlot(corTab,params,flags)  
    figure; hold on;
    uAges = unique(corTab.age);
    tmpTab = table();
    tmpTab.ranges = unique(corTab.range);
    axs = [];
    for ii = 1:length(uAges)
        ii;
        age = uAges(ii);
        subTab = corTab(corTab.age(:,1) == age,:);
        tmpMeans = []; tmpStd = [];
        for jj = 1:height(tmpTab)
            tmpMeans = [tmpMeans nanmean(subTab.cors(subTab.range == tmpTab.ranges(jj)))];
            tmpStd = [tmpStd nanstd(subTab.cors(subTab.range == tmpTab.ranges(jj)))/sqrt(length(subTab.cors(subTab.range == tmpTab.ranges(jj))))];
        end
        h = errorbar(tmpTab.ranges+0.05*ii - 0.25,tmpMeans,tmpStd,'LineWidth',2);
        axs = [axs h];
    end
    poolMean = []; poolStd = [];
    for jj = 1:height(tmpTab)
        poolMean = [poolMean nanmean(corTab.cors(corTab.range == jj))];
        poolStd = [poolStd nanstd(corTab.cors(corTab.range == jj))/sqrt(length(corTab.cors(corTab.range == jj)))];
        for kk = jj:height(tmpTab)
            if flags.meanCorr %length(corTab.cors(corTab.range == jj)) == length(corTab.cors(corTab.range == kk))
                [h,p,ci,stats] = ttest(corTab.cors(corTab.range == kk),corTab.cors(corTab.range == jj));
                if h == 1        
                   sigstar({[kk,jj]},p*factorial(height(tmpTab)))
                end
            else
                [h,p,ci,stats] = ttest2(corTab.cors(corTab.range == kk),corTab.cors(corTab.range == jj));
                if h == 1        
                   sigstar({[kk,jj]},p*factorial(height(tmpTab)))
                end                
            end
        end
    end
    h = errorbar(tmpTab.ranges+0.05*ii - 0.2,poolMean,poolStd,'LineWidth',2,'color',rgb('black'));
    axs = [axs h];
    xticks([1:4])
    xticklabels({sprintf('Corr %d-%d \\mu m', params.lNeigh , params.mNeigh),...
        sprintf('Corr %d-%d \\mu m', params.mNeigh , params.umNeigh),...
        sprintf('Corr %d-%d \\mu m', params.umNeigh , params.uNeigh),...
        sprintf('Corr > %d \\mu m', params.uNeigh)})
    ylabel('correlation coefficient')
    xlim([0.5,4.5]); ylim([0,1])
    legend(axs,{'P8','P9','P10','P11','P14','pooled'},'Location','southwest')
    %set(gca,'FontSize',25) 
end

function [] = cpdPlot(cpdTab)
    figure;
    dx = linspace(0,60,300);
    subplot(2,4,1)
    hold on
    plot(dx,cpdTab.n1(1:300)/sum(cpdTab.n1(1:300)),'LineWidth',2)
    plot(dx,cpdTab.n1(1:300)/sum(cpdTab.n1(1:300)) + sqrt(cpdTab.n1(1:300))/sum(cpdTab.n1(1:300)),'r--','LineWidth',1)
    plot(dx,cpdTab.n1(1:300)/sum(cpdTab.n1(1:300)) - sqrt(cpdTab.n1(1:300))/sum(cpdTab.n1(1:300)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])

    title('5-75 nm')
    ylabel('Probability')
    subplot(2,4,2)
    hold on
    plot(dx,cpdTab.n2(1:300)/sum(cpdTab.n2(1:300)),'LineWidth',2)
    plot(dx,cpdTab.n2(1:300)/sum(cpdTab.n2(1:300)) + sqrt(cpdTab.n2(1:300))/sum(cpdTab.n2(1:300)),'r--','LineWidth',1)
    plot(dx,cpdTab.n2(1:300)/sum(cpdTab.n2(1:300)) - sqrt(cpdTab.n2(1:300))/sum(cpdTab.n2(1:300)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('75-150 nm')
    xlabel('Time (s)')
    subplot(2,4,3)
    hold on
    plot(dx,cpdTab.n3(1:300)/sum(cpdTab.n3(1:300)),'LineWidth',2)
    plot(dx,cpdTab.n3(1:300)/sum(cpdTab.n3(1:300)) + sqrt(cpdTab.n3(1:300))/sum(cpdTab.n3(1:300)),'r--','LineWidth',1)
    plot(dx,cpdTab.n3(1:300)/sum(cpdTab.n3(1:300)) - sqrt(cpdTab.n3(1:300))/sum(cpdTab.n3(1:300)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('150-225 nm')
    subplot(2,4,4)
    hold on
    plot(dx,cpdTab.n4(1:300)/sum(cpdTab.n4(1:300)),'LineWidth',2)
    plot(dx,cpdTab.n4(1:300)/sum(cpdTab.n4(1:300)) + sqrt(cpdTab.n4(1:300))/sum(cpdTab.n4(1:300)),'r--','LineWidth',1)
    plot(dx,cpdTab.n4(1:300)/sum(cpdTab.n4(1:300)) - sqrt(cpdTab.n4(1:300))/sum(cpdTab.n4(1:300)),'r--','LineWidth',1)
    xlim([0,15]);
    title('>225 nm')
    yyaxis right
    ylabel('L events')
    yticks([])


    dx = linspace(0,60,300);
    subplot(2,4,5)
    hold on
    plot(dx,cpdTab.n1(301:end)/sum(cpdTab.n1(301:end-2)),'LineWidth',2)
    plot(dx,cpdTab.n1(301:end)/sum(cpdTab.n1(301:end-2)) + sqrt(cpdTab.n1(301:end))/sum(cpdTab.n1(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n1(301:end)/sum(cpdTab.n1(301:end-2)) - sqrt(cpdTab.n1(301:end))/sum(cpdTab.n1(301:end-2)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])

    title('5-75 nm')
    ylabel('Probability')
    subplot(2,4,6)
    hold on
    plot(dx,cpdTab.n2(301:end)/sum(cpdTab.n2(301:end-2)),'LineWidth',2)
    plot(dx,cpdTab.n2(301:end)/sum(cpdTab.n2(301:end-2)) + sqrt(cpdTab.n2(301:end))/sum(cpdTab.n2(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n2(301:end)/sum(cpdTab.n2(301:end-2)) - sqrt(cpdTab.n2(301:end))/sum(cpdTab.n2(301:end-2)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('75-150 nm')
    xlabel('Time (s)')
    subplot(2,4,7)
    hold on
    plot(dx,cpdTab.n3(301:end)/sum(cpdTab.n3(301:end-2)),'LineWidth',2)
    plot(dx,cpdTab.n3(301:end)/sum(cpdTab.n3(301:end-2)) + sqrt(cpdTab.n3(301:end))/sum(cpdTab.n3(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n3(301:end)/sum(cpdTab.n3(301:end-2)) - sqrt(cpdTab.n3(301:end))/sum(cpdTab.n3(301:end-2)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('150-225 nm')
    subplot(2,4,8)
    hold on
    plot(dx,cpdTab.n4(301:end)/sum(cpdTab.n4(301:end)),'LineWidth',2)
    plot(dx,cpdTab.n4(301:end)/sum(cpdTab.n4(301:end-2)) + sqrt(cpdTab.n4(301:end))/sum(cpdTab.n4(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n4(301:end)/sum(cpdTab.n4(301:end-2)) - sqrt(cpdTab.n4(301:end))/sum(cpdTab.n4(301:end-2)),'r--','LineWidth',1)
    xlim([0,15]);
    title('>225 nm')
    yyaxis right
    ylabel('H events')
    yticks([])


    %%
    savTab = cpdTab
    %%
    cpdTab(1:end-2,:) = array2table(table2array(savTab(1:end-2,:)) );
    %%
    prH = table2array(cpdTab(end,1:3))/sum(table2array(cpdTab(end,1:3)));
    prL = table2array(cpdTab(end-1,1:3))/sum(table2array(cpdTab(end-1,1:3)));


    Htab = table2array(cpdTab(1:300,1:3)) ;
    Ltab = table2array(cpdTab(301:end-2,1:3));

    pdt_rH = Htab./repmat(sum(Htab,1),[300,1]);
    pdt_rL = Ltab./repmat(sum(Ltab,1),[298,1]);
    pdtH = pdt_rH*prH';
    pdtL = pdt_rL*prL';

    IRH = pdt_rH.*log2((pdt_rH )./repmat(pdtH , [1,3]));
    IRH(isnan(IRH)) = 0;
    IRH = sum(IRH,1)*prH';

    IRH_corrected = IRH - (size(Htab,1)*size(Htab,2)-size(Htab,1)-size(Htab,2))/(2*sum(table2array(cpdTab(end,1:3)))*log(2));

    IRL = pdt_rL.*log2((pdt_rL)./repmat(pdtL , [1,3]));
    IRL(isnan(IRL)) = 0;
    IRL = sum(IRL,1)*prL';
    IRL_corrected = IRL - (size(Htab,1)*size(Htab,2)-size(Htab,1)-size(Htab,2))/(2*sum(table2array(cpdTab(end-1,1:3)))*log(2));


    %%

    pcpdTab = table2array(cpdTab(1:300,:)) + table2array(cpdTab(301:end,:));

    pcpdTab(299:300,:) = 0;

    figure;
    dx = linspace(0,60,300);
    subplot(1,3,1)
    hold on
    plot(dx,pcpdTab(:,1)/sum(pcpdTab(:,1)),'LineWidth',2)
    ylabel('Probability')
    title('5-35 nm')
    ylim([0,0.13])
    subplot(1,3,2)
    hold on
    plot(dx,pcpdTab(:,2)/sum(pcpdTab(:,2)),'LineWidth',2)
    xlabel('Time (s)')
    title('35-105 nm')
    ylim([0,0.13])

    subplot(1,3,3)
    hold on
    plot(dx,pcpdTab(:,3)/sum(pcpdTab(:,3)),'LineWidth',2)
    title('105-245 nm')
    ylim([0,0.13])

    %%

    xlim([0,15]);ylim([0,0.5])

    title('5-75 nm')
    ylabel('Probability')
    subplot(2,4,2)
    hold on
    plot(dx,cpdTab.n2(1:300)/sum(cpdTab.n2(1:300)),'LineWidth',2)
    plot(dx,cpdTab.n2(1:300)/sum(cpdTab.n2(1:300)) + sqrt(cpdTab.n2(1:300))/sum(cpdTab.n2(1:300)),'r--','LineWidth',1)
    plot(dx,cpdTab.n2(1:300)/sum(cpdTab.n2(1:300)) - sqrt(cpdTab.n2(1:300))/sum(cpdTab.n2(1:300)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('75-150 nm')
    xlabel('Time (s)')
    subplot(2,4,3)
    hold on
    plot(dx,cpdTab.n3(1:300)/sum(cpdTab.n3(1:300)),'LineWidth',2)
    plot(dx,cpdTab.n3(1:300)/sum(cpdTab.n3(1:300)) + sqrt(cpdTab.n3(1:300))/sum(cpdTab.n3(1:300)),'r--','LineWidth',1)
    plot(dx,cpdTab.n3(1:300)/sum(cpdTab.n3(1:300)) - sqrt(cpdTab.n3(1:300))/sum(cpdTab.n3(1:300)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('150-225 nm')
    subplot(2,4,4)
    hold on
    plot(dx,cpdTab.n4(1:300)/sum(cpdTab.n4(1:300)),'LineWidth',2)
    plot(dx,cpdTab.n4(1:300)/sum(cpdTab.n4(1:300)) + sqrt(cpdTab.n4(1:300))/sum(cpdTab.n4(1:300)),'r--','LineWidth',1)
    plot(dx,cpdTab.n4(1:300)/sum(cpdTab.n4(1:300)) - sqrt(cpdTab.n4(1:300))/sum(cpdTab.n4(1:300)),'r--','LineWidth',1)
    xlim([0,15]);
    title('>225 nm')
    yyaxis right
    ylabel('L events')
    yticks([])


    dx = linspace(0,60,300);
    subplot(2,4,5)
    hold on
    plot(dx,cpdTab.n1(301:end)/sum(cpdTab.n1(301:end-2)),'LineWidth',2)
    plot(dx,cpdTab.n1(301:end)/sum(cpdTab.n1(301:end-2)) + sqrt(cpdTab.n1(301:end))/sum(cpdTab.n1(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n1(301:end)/sum(cpdTab.n1(301:end-2)) - sqrt(cpdTab.n1(301:end))/sum(cpdTab.n1(301:end-2)),'r--','LineWidth',1)
    %xlim([0,15]);ylim([0,0.5])

    title('5-75 nm')
    ylabel('Probability')
    subplot(2,4,6)
    hold on
    plot(dx,cpdTab.n2(301:end)/sum(cpdTab.n2(301:end-2)),'LineWidth',2)
    plot(dx,cpdTab.n2(301:end)/sum(cpdTab.n2(301:end-2)) + sqrt(cpdTab.n2(301:end))/sum(cpdTab.n2(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n2(301:end)/sum(cpdTab.n2(301:end-2)) - sqrt(cpdTab.n2(301:end))/sum(cpdTab.n2(301:end-2)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('75-150 nm')
    xlabel('Time (s)')
    subplot(2,4,7)
    hold on
    plot(dx,cpdTab.n3(301:end)/sum(cpdTab.n3(301:end-2)),'LineWidth',2)
    plot(dx,cpdTab.n3(301:end)/sum(cpdTab.n3(301:end-2)) + sqrt(cpdTab.n3(301:end))/sum(cpdTab.n3(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n3(301:end)/sum(cpdTab.n3(301:end-2)) - sqrt(cpdTab.n3(301:end))/sum(cpdTab.n3(301:end-2)),'r--','LineWidth',1)
    xlim([0,15]);ylim([0,0.5])
    title('150-225 nm')
    subplot(2,4,8)
    hold on
    plot(dx,cpdTab.n4(301:end)/sum(cpdTab.n4(301:end)),'LineWidth',2)
    plot(dx,cpdTab.n4(301:end)/sum(cpdTab.n4(301:end-2)) + sqrt(cpdTab.n4(301:end))/sum(cpdTab.n4(301:end-2)),'r--','LineWidth',1)
    plot(dx,cpdTab.n4(301:end)/sum(cpdTab.n4(301:end-2)) - sqrt(cpdTab.n4(301:end))/sum(cpdTab.n4(301:end-2)),'r--','LineWidth',1)
    xlim([0,15]);
    title('>225 nm')
    yyaxis right
    ylabel('H events')
    yticks([])

end