%%
% Combines the results from several trials of a single experiment. All
% .csv's must be collected in one folder.
%   Input: pathTo - path to the folder containing the csv's
%          fixThresh - set this value if you want to fix the threshold to
%                      the quantile = fixThresh
%   Output: tab - a table containing:
%                       -event rates
%                       -event durations
%                       -event start times
%                       -event amplitudes
%                       -event jitter
%

function [tab,pcaTab,corTab,cpdTab] = combineFolder(pathTo, params , flags)
      
    fixThresh = params.fixThresh;
    % Siegel proposes to cutoff at 20%
    LOWERBOUND = params.LOWERBOUND;
    % get all files in the folder
    fileList = getFileList(pathTo);
    
    % initialize arrays
    tab = array2table(nan(1,length(params.outTableNames)),'VariableNames',params.outTableNames); 
    pcaTab = table(); corTab = table();
    cpdTab = zeros(600,4);
    
    cID = NaN; cCond = NaN;
    
    LHrat = [];
    
    for ii = 1:length(fileList)
        fileStr = fileList{ii}
        % extract expriment identifiers
        [expID , expID2] = getExpID(fileStr,params);

        % get the average activity from the ROI's.
        [readIn ,skip] = importCSV([pathTo fileStr] , params , flags);
        % If bad animal or too many ROIs
        if skip || size(readIn.means,2) > params.upperN
            continue 
        end
        % extract thresholds
        if flags.fixThreshPerExp
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
        else
            [thrTraces , thrDiffTraces] = get_Threshold(readIn.means,fixThresh,flags); 
        end
        % convert to raster plot
        [raster , ~] = TracesToSpikeTimes(readIn.means , thrTraces , params.trDiffRat*thrDiffTraces);
        T = readIn.dT * size(raster,1);
        % PCA analysis for 80% of variance
        [pcaPerc , pcaPercAlt ] = getPcaPerc(raster , params);

        % get naive participation rate
        [ParRate , duds] = get_ParRate(raster , params.nParRateWindow , params.dudslim);
        
        % convert to event rate 0.2
        [evRates , evDurations ,evDurationsFWHM , evStartTime] = get_evRates(raster , ParRate ,...
            LOWERBOUND , prctile(diff(ParRate(:)),params.diffEvRatePCT) , ceil(params.evRateWindow*readIn.dT/0.13) , duds , params);
        if flags.plotRast
            plotRast(raster,readIn.means,evStartTime,evDurations,evRates,ParRate , fileStr , params,flags);
        end%smoothdata(raster,1,'gaussian',params.pcaSmoothWindowSize),
        
        % extract amplitudes and jitter
        if flags.ampF0Flag
            evAmps = get_evAmps(readIn.means , raster , evStartTime , evDurations);%evDurations);        
        else
            evAmps = get_evAmps(readIn.sMeans , raster , evStartTime , evDurations);%evDurations);
        end
        evJitter = get_evJitter(raster , evStartTime , evDurations).*(evDurationsFWHM./evDurations);

        if flags.plotLHpos
            
            figure(); hold on;
            scatter(readIn.pos(:,2),-readIn.pos(:,1),50,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
            axis equal
            LHratColl = zeros(size(readIn.pos,1) , length(evStartTime));
            for jj = 1:length(evStartTime)
                if jj > 1 & evRates(jj-1) > 0.8
                    figure(); hold on;
                    scatter(readIn.pos(:,2),-readIn.pos(:,1),50,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
                    axis equal
                end

                if evRates(jj) > 0.8
                    parti = find(sum(raster(evStartTime(jj):evStartTime(jj)+evDurations(jj) , :),1) > 0 );
                    partiA = max(readIn.means(evStartTime(jj):evStartTime(jj)+evDurations(jj) , parti))./max(readIn.means(: , parti));
                    LHratColl(parti , jj) = partiA;
                    scatter(readIn.pos(parti,2) + (rand(length(parti),1) -0.5)*5 ,...
                        -readIn.pos(parti,1) + (rand(length(parti),1) -0.5)*5,partiA.^2*70,'Marker','*','MarkerFaceColor','blue','MarkerEdgeColor','blue')
                    LHrat = [LHrat ; partiA./mean(LHratColl(parti , 1:jj-1),2)];
                else
                    parti = find(sum(raster(evStartTime(jj):evStartTime(jj)+evDurations(jj) , :),1) > 0 );
                    partiB = max(readIn.means(evStartTime(jj):evStartTime(jj)+evDurations(jj) , parti))./max(readIn.means(: , parti));
                    LHratColl(parti , jj) = partiB;
                    scatter(readIn.pos(parti,2) + (rand(length(parti),1) -0.5)*5 ,...
                        -readIn.pos(parti,1) + (rand(length(parti),1) -0.5)*5,partiB.^2*70,'Marker','o','MarkerFaceAlpha',0.5,'MarkerFaceColor','red','MarkerEdgeColor','black')
                end
            end
            
        end
        
        if flags.corrFlag || flags.cpdFlag
            sRaster = smoothdata(raster,1,'gaussian',params.pcaSmoothWindowSize);
            [adjMat , corMat] = getAdjCor(sRaster, readIn.pos , readIn.dimReal , readIn.dimPix);
            if flags.storeCorr
                [~,fName,~] = fileparts(fileStr);
                save(['mat/' fName(params.plotRastl:end) '.mat'] , 'sRaster' , 'readIn');
            end
            if  ~isempty(adjMat) && flags.cpdFlag
                cpdTab = cpdTab + getCpd(raster , evStartTime , evDurations , evRates , adjMat , readIn.dT ,  params);                
            end
            if ~isempty(adjMat) && ~isempty(corMat) && flags.corrFlag
                if flags.corrPlot
                    f = figure('units','normalized','outerposition',[0 0 1 1],'Visible',flags.plotVisible);
                    subplot(1,2,1),
                    imshow(adjMat , [])
                    subplot(1,2,2)
                    imshow(corMat , [])
                    [~,fName,~] = fileparts(fileStr);
                    print(['fig_' params.Dataset '/corMat/' fName(params.ID1s:end) '.png'],'-dpng')
                    close(f);
                end
                corTab = [corTab ; getCorTab(adjMat, corMat , readIn.age , expID  , readIn.cond ,  params,flags)];
            end
        end 
        
        
        tab = [tab ; table(evRates ,evDurations*readIn.dT ,evDurationsFWHM*readIn.dT ,...
                    transp(evStartTime) ,evAmps ,transp(evStartTime).*readIn.dT,...
                    evJitter*readIn.dT , repmat(readIn.cond , size(evJitter)) , ...
                    repmat(readIn.F0 , size(evJitter)) , repmat(readIn.F0quant , size(evJitter)) , ...
                    repmat(readIn.age , size(evJitter)) , repmat(expID , size(evJitter)) , ...
                    repmat(expID2 , size(evJitter)) , repmat(T , size(evJitter)) ,...
                    'VariableNames',params.outTableNames)];
        tab = [tab ; array2table(nan(1,length(params.outTableNames)),'VariableNames',params.outTableNames)];
        
        pcaTab = [pcaTab; table( pcaPerc' , pcaPercAlt' , repmat(readIn.age ,length(pcaPerc),1) ,...
            repmat(expID,length(pcaPerc),1) , repmat(expID2,length(pcaPerc),1),...
            repmat(readIn.cond,length(pcaPerc),1) , 'VariableNames',params.pcaTableNames)];

    end
end

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

function [pcaPerc , pcaPercAlt ]= getPcaPerc(raster , params)
        tempSmoothRaster = smoothdata(raster,1,'gaussian',params.pcaSmoothWindowSize);%raster'*(eye(size(raster,1),size(raster,1)) +diag(ones(size(raster,1)-1,1),1) + diag(ones(size(raster,1)-2,1),2));
        [~ , ~ , ~ , ~, exp , ~] = pca(tempSmoothRaster);
        
        cumExp = cumsum(exp);
        [~,idx] = max(cumExp > params.PCAbound);
        pcaPerc = idx*1.0/length(exp);
        pcaPercAlt = size(raster,2)*1.0/sum(exp.^2);
end

function [] = plotRast(raster,roiVals,evStartTime,evDurations,evRates,ParRate , file , params,flags)
    f = figure('units','normalized','outerposition',[0 0 1 1],'visible',flags.plotVisible);
    subtightplot(2,1,1)
    %plot(smoothdata(roiVals)./repmat(max(roiVals,[],1),size(roiVals,1),1) +repmat(size(roiVals,2):-1:1,size(roiVals,1),1),'-','color',rgb('black'))
    plot(roiVals./repmat(max(roiVals,[],1),size(roiVals,1),1) +repmat(size(roiVals,2):-1:1,size(roiVals,1),1),'-','color',rgb('black'))

    xlim([0, size(roiVals,1)]); ylim([0 , size(roiVals,2)+1])
    subtightplot(2,1,2)

    imagesc(raster','CDataMapping','scaled');
    hold on; plot(ParRate*size(raster,2));
    %scatter(evStartTime , ones(length(evStartTime),1)*size(raster,2)*0.5,'filled');
    scatter(evStartTime , evRates*size(raster,2) , 'filled')
    scatter(evStartTime+evDurations' , evRates*size(raster,2) , 'filled')

    [~,fName,~] = fileparts(file);
    colormap(jet)
    % print(['fig_' params.Dataset '/tracesAndRaster/' fName(params.ID1s:end) '.svg'],'-dsvg')
    close(f)

%     f = figure;
%     RasterWrapper(raster,evStartTime,evDurations,evRates,ParRate);
%     [~,fName,~] = fileparts(file);
%     print(['fig_Fried/rasters/' fName(params.plotRastl:end) '.png'],'-dpng')
%     close(f)
end

function  [adjMat , corMat] = getAdjCor(roiVals , pos , dimReal , dimPix)
            dimRatio = dimReal./dimPix;
            microPos = pos.*repmat(dimRatio,size(pos,1),1);

            adjMat = getAdjMat(microPos);

            corMat = corr(roiVals);
end

function  cpdTab = getCpd(raster , evStartTime , evDurations , evRates , adjMat , dT ,  params)
    lowerNeighbor = params.lNeigh; midNeighbor = params.mNeigh;
    upperMidNeighbor = params.umNeigh; upperneighbor = params.uNeigh;

    nearCPD = zeros(600,1);
    midCPD = zeros(600,1);
    uppermidCPD = zeros(600,1);
    farCPD = zeros(600,1);
    for kk1 = 1:size(raster,2)
        for kk2p = 1:length(evStartTime)%size(raster,1)
            kk2 = evStartTime(kk2p);
            pp2 = evRates(kk2p);
            for dur = 0:evDurations(kk2p)
            if raster(max([1,kk2+dur-1]),kk1)
                for kk3 = 1:size(raster,2)
                    if adjMat(kk1,kk3) > lowerNeighbor && adjMat(kk1,kk3) < midNeighbor && kk1~= kk3
                        [maxVal,maxInd] = max(raster(kk2:min(kk2+99,size(raster,1)),kk3));
                        if maxVal*maxInd > 0
                            nearCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) = nearCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) +1;
                            if pp2 > 0.8
                                nearCPD(end) = nearCPD(end) +1;
                            else
                                nearCPD(end-1) = nearCPD(end-1) +1;                                            
                            end
                        end
                    elseif adjMat(kk1,kk3) > midNeighbor && adjMat(kk1,kk3) < upperMidNeighbor && kk1~= kk3
                        [maxVal,maxInd] = max(raster(kk2:min(kk2+99,size(raster,1)),kk3));
                        if maxVal*maxInd > 0
                            midCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) = midCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) +1;
                            if pp2 > 0.8
                                midCPD(end) = midCPD(end) +1;
                            else
                                midCPD(end-1) = midCPD(end-1) +1;                                            
                            end
                        end
                    elseif adjMat(kk1,kk3) > upperMidNeighbor && adjMat(kk1,kk3) < upperneighbor && kk1~= kk3
                        [maxVal,maxInd] = max(raster(kk2:min(kk2+99,size(raster,1)),kk3));
                        if maxVal*maxInd > 0
                            uppermidCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) = uppermidCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) +1;
                            if pp2 > 0.8
                                uppermidCPD(end) = uppermidCPD(end) +1;
                            else
                                uppermidCPD(end-1) = uppermidCPD(end-1) +1;
                           end
                        end
                    elseif adjMat(kk1,kk3) > upperneighbor && kk1~= kk3
                        [maxVal,maxInd] = max(raster(kk2:min(kk2+99,size(raster,1)),kk3));
                        if maxVal*maxInd > 0
                            farCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) = farCPD(ceil(maxInd*dT*5) + 300*(pp2 > 0.8)) +1;
                            if pp2 > 0.8
                                farCPD(end) = farCPD(end) +1;                                    
                            else
                                farCPD(end-1) = farCPD(end-1) +1;                                    
                            end
                        end
                    end
                end
            end
            end
        end
    end  
    cpdTab = [nearCPD ,midCPD,uppermidCPD,farCPD];
end

function corTab = getCorTab(adjMat, corMat , age , expID , cond , params,flags)
    lowerNeighbor = params.lNeigh; midNeighbor = params.mNeigh;
    upperMidNeighbor = params.umNeigh; upperneighbor = params.uNeigh;
    % Correlation analysis
    neighMask = (adjMat > lowerNeighbor) .* (adjMat < midNeighbor);
    neighMask(neighMask == 0) = NaN;
    neighCor = corMat.*neighMask;
    neighCor1 = neighCor(~isnan(neighCor));

    neighMask = (adjMat > midNeighbor) .* (adjMat < upperMidNeighbor);
    neighMask(neighMask == 0) = NaN;
    neighCor = corMat.*neighMask;
    neighCor2 = neighCor(~isnan(neighCor));

    neighMask = (adjMat > upperMidNeighbor) .* (adjMat < upperneighbor);
    neighMask(neighMask == 0) = NaN;
    neighCor = corMat.*neighMask;
    neighCor3 = neighCor(~isnan(neighCor));

    neighMask = (adjMat > upperneighbor) *1.0;
    neighMask(neighMask == 0) = NaN;
    neighCor = corMat.*neighMask;
    neighCor4 = neighCor(~isnan(neighCor));

    corTabTmp = table();
    if flags.meanCorr
        corTabTmp.cors = [mean(neighCor1) ; mean(neighCor2) ; mean(neighCor3) ; mean(neighCor4)];
        corTabTmp.range = [1;2;3;4];%ones(length(neighCor1),1);2*ones(length(neighCor2),1);3*ones(length(neighCor3),1);4*ones(length(neighCor4),1)];
        corTabTmp.age = repmat(age,4,1);%length([neighCor1 ; neighCor2 ; neighCor3 ; neighCor4]),1);
        corTabTmp.id = repmat(expID,4,1);%length(corTabTmp.age),1);
        corTabTmp.cond = repmat(cond,4,1);%length(corTabTmp.age),1);

    else
        corTabTmp.cors = [neighCor1 ; neighCor2 ; neighCor3 ; neighCor4];
        corTabTmp.range = [ones(length(neighCor1),1);2*ones(length(neighCor2),1);3*ones(length(neighCor3),1);4*ones(length(neighCor4),1)];
        corTabTmp.age = repmat(age,length([neighCor1 ; neighCor2 ; neighCor3 ; neighCor4]),1);
        corTabTmp.id = repmat(expID,length(corTabTmp.age),1);
        corTabTmp.cond = repmat(cond,length(corTabTmp.age),1);

    end
    corTab = corTabTmp;% [corTab ; corTabTmp];
end