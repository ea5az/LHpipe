%%%
% Imports a .csv file created from the Multi-Measure tool of ImageJ/Fiji.
% Input : csvPath - Path to the .csv file
%         correctF0 - compute F/F0 flag
%         filterOut - flag for filtering out bad ROIs
function [readIn ,skip] = importCSV(csvPath , params , flags)

    % Extract flags
    if nargin < 3
        correctF0 = 0; filterOut = 0;
    else
        correctF0 = flags.correctF0;
        filterOut = flags.filterOut;
    end
    
    % Get persitent experiment table
    persistent fileTab;
    if isempty(fileTab)
       tabPath = params.tabPath;
       fileTab =  readtable(tabPath);
    end
    
    % Initialize arrays
    ROIpath = params.ROIpath;
    readIn = struct;
    readIn.means = [];readIn.maxs = [];readIn.mins = [];readIn.dT = 0;
    readIn.dimPix = [];readIn.dimReal = [];readIn.T = 0; 
    readIn.N = 0; readIn.cond = ''; readIn.age = 0;readIn.skip = 0;
    readIn.F0=NaN;readIn.F0quant = 0;readIn.pos=[];
    
    % Extract experiment identifiers
    [idx , expName] = getIDX(fileTab , csvPath , params);
    
    %Skip if animal not found or subprime animals are excluded
    if isempty(idx)
        skip = 1; return;
    end
       
    % Get experiment condition from table
    [cond , skip] = getCond(fileTab , idx , params , flags);
    
    if flags.onlyPrime
        if ~strcmp(fileTab.ISPRIME{idx} , 'PRIME')
            skip = 1; return;
        end
    end

    if skip
        return;
    end
    
    % Set condition and ROI positions
    readIn.cond = cond;
    readIn.pos = readPos(ROIpath , expName);
    
    % scatter ROI positions
    if flags.scatterPos
        f = figure(); hold on
        scatter(readIn.pos(:,2),-readIn.pos(:,1))
    end
    
    % If number of cells provided in experiment sheet
    if ~isempty(fileTab.x_Cells{idx}) && ~strcmp(fileTab.x_Cells{idx},'?')
        readIn.N = str2num(fileTab.x_Cells{idx});
    else % otherwise assume that experiment was subprime
        skip = 1; 
        if flags.scatterPos; close(f); end
        return;
    end
    
    % import traces from file
    [T , means , ~ , ~] = getTraces(csvPath);
    readIn.sMeans = means;
    
    % check table if portion of the experiment should be excluded
    tlim = fileTab.length{idx};
    if ~isempty(tlim) && ~strcmp(tlim,'?')
        tlim = strsplit(tlim,{':','-'});
        t1 = str2double(tlim{1});
        t2 = str2double(tlim{2});
        means = means(t1:t2,:);
        T = t2-t1;
        readIn.sMeans = means;
    end
    % extract recording FOV dimension
    [dimPix , dimReal] = getDim(fileTab , idx);
    % exclude ROIS at the edge of FOV
    if filterOut
        edgeMask = filterEdges(readIn.pos , dimPix , params);
        F0means = means(:,edgeMask); readIn.F0quant = mean(F0means(:));
        means(:,edgeMask) = []; readIn.pos(edgeMask,:) = [];
        mask = filterMeans(means , params , readIn.N );
        % leave out certain percentage of ROIs
        if sum(mask) > 0
            F0means = means(:,mask); readIn.F0quant = mean(F0means(:));
            means(:,mask) = []; readIn.pos(mask,:) = [];
        end
        scatter(readIn.pos(:,2),-readIn.pos(:,1),'filled')
        idMat = pdist2(readIn.pos,readIn.pos)+eye(size(readIn.pos,1))*100 < params.idLim;
        idMat = idMat - triu(idMat);
        [row , col] = find(idMat);
        cMeans = corr(means,means);
        cThresh = prctile(cMeans(:),params.corrThresh);
        if length(row) > 1
            for ii = 1:length(row)
                corr(means(:,row(ii)),means(:,col(ii)))
                if corr(means(:,row(ii)),means(:,col(ii))) > cThresh
                    means(:,row(ii)) = []; readIn.pos(row(ii),:) = [];
                end
            end
        else
            corr(means(:,row),means(:,col))
            if corr(means(:,row),means(:,col)) > cThresh
                means(:,row) = []; readIn.pos(row,:) = []; 
            end
        end        

        readIn.sMeans = means;
    end
    % add reduced ROIs in different color
    if flags.scatterPos
        scatter(readIn.pos(:,2),-readIn.pos(:,1),'filled')
        axis equal
%         print(['fig_Fried/ROI/' expName '.png'],'-dpng')
        close(f)
    end
    % normalize each cell by its own F0
    if correctF0 
        a = params.FilterA;
        means = (means - repmat(mean(means,1),size(means,1),1))./repmat(mean(means,1),size(means,1),1);%mean(F0means(:));
        readIn.sMeans = means;
        % low pass filter
        if flags.lowFilt
            means = filter(a, [1 a-1], means);
        end
        % high pass filter
        if flags.highFilt
            means = filter([1-a/2 a/2-1],[1 a/2-1], means);
        end
    end
    
    if flags.detrend
       means =  detrend(means);
    end
    % extract age
    ageStr = fileTab.age{idx};
    if ~strcmp(ageStr,'?')
        age = str2num(ageStr(2:end));
    else
        age = 0; skip = 1; return;
    end
    % store everything in structure
    readIn.F0 = fileTab.F0(idx); readIn.T = T;
    readIn.fr = fileTab.frameRate{idx}; readIn.dT = str2double(readIn.fr(1:end-1));
    readIn.means = means; 
    readIn.age = age; readIn.dimPix = dimPix; readIn.dimReal = dimReal;
end

function edgeMask = filterEdges(pos , dimPix , params)
    edgeMask =(pos(:,1) <= params.EdgePix) + (pos(:,2) <= params.EdgePix) + ...
        (abs(pos(:,1) - dimPix(2)) <= params.EdgePix) + (abs(pos(:,2) - dimPix(1)) <= params.EdgePix);
    edgeMask = edgeMask > 0;
end

function pos = readPos(ROIpath , expName)
    posTab = readtable(fullfile(ROIpath,['ROI_MOCO_' expName '.csv']));
    pos = table2array(posTab);
end

function [idx , expName] = getIDX(fileTab , csvPath , params)
    [~,fullExpName,~] = fileparts(csvPath);
    expName = fullExpName(params.ExpNameStart:end);
    idx =find(strcmp(fileTab.Var1,expName));
end

function [cond , skip] = getCond(fileTab , idx , params , flags)
    condStr = fileTab.Condition{idx};
    conditions = params.conditions;
    cond = length(conditions);
    skip = 0;
    for ii = 1:length(conditions)
        if strcmp(condStr , conditions{ii})
            cond = ii - 1;
            ssplit = strsplit(conditions{ii});
            ssplit = ssplit{1};
            if flags.(['skip' ssplit])
                skip = 1; 
            end
        end
    end
    if cond == length(conditions)
        skip = 1;
    end
end

function [T , means , maxs , mins] = getTraces(csvPath)
    tab = readtable(csvPath,'TreatAsEmpty',{'Infinity','-Infinity'});
    colNames = tab.Properties.VariableNames;
    
    T = tab.Var1;
    colNames = colNames(2:end);
    
    means = zeros(length(T) , length(colNames)/4);
    maxs = zeros(length(T) , length(colNames)/4);
    mins = zeros(length(T) , length(colNames)/4);
    for ii=1:length(colNames)/4
        means(:,ii) = table2array(tab(:,4*(ii-1) + 3));
        maxs(:,ii) = table2array(tab(:,4*(ii-1) + 4));
        mins(:,ii) = table2array(tab(:,4*(ii-1) + 5));
    end
end

function mask = filterMeans(means , params , N)
    varMeans = var(means);
    sortVar = sort(varMeans,'descend');
    bound = sortVar(min(length(sortVar),floor(params.leaveOutRatio*N)+1));
    mask =  varMeans < bound;
end

function [dimPix , dimReal] = getDim(fileTab , idx)
    dimPix = cell2mat(strsplit(fileTab.FOV_pix_{idx},'x')');
    dimPix = [str2num(dimPix(1,:)) str2num(dimPix(2,:))];
    dimReal = cell2mat(strsplit(fileTab.FOV_um_{idx},'x')');
    dimReal = [str2num(dimReal(1,:)) str2num(dimReal(2,:))];
end