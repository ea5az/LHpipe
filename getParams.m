function params = getParams(DATASET)
    params = struct();
    params.Dataset = DATASET;
    if strcmp(DATASET , 'Fried')
        % Set paths
        params.ROIpath = '/storage/gjor/Data/kirchnerj/outDirFried/mat/';
        params.tabPath = '/storage/gjor/Data/kirchnerj/Data/Friederike.csv';
        params.TRACESpath = '/storage/gjor/Data/kirchnerj/Data/Friederike/TRACES/';
        %params.TRACESpath = '/storage/gjor/Data/kirchnerj/Data/Friederike/TEST/';
        % conditions in this study
        params.conditions = {'NORMAL' , 'BINOCULAR ENUCLEATION' , 'FORSKOLIN INJECTION (BINOCULAR)' , 'SUBPAR'};
        % formating strings
        params.ExpNameStart = 12;
        params.plotRastl = 12;
        params.ID1s = 16;
        params.ID1e = 7;
        params.ID2s = 19;
        params.ID2e = 4;
        
        params.fixThresh = 94;%96; % percentile of voltage threshold
        params.diffEvRatePCT = 97; % percentile of participation threshold
        params.leaveOutRatio = 1.;%0.95;%0.9; % how many percent of ROIs to keep
        params.trDiffRat = 1;
        params.endEventPR = 1.7; % multiply event threshold by this number to get end threshold
        params.evRateWindow = 21; % how many events to consider maximally

        params.labels = {'NORMAL' , 'BINOCULAR ENUCLEATION' , 'FORSKOLIN INJECTION (BINOCULAR)' };

    elseif strcmp(DATASET , 'Curly')
        params.ROIpath = '/storage/gjor/Data/kirchnerj/outDirCurly/mat/';
        params.tabPath = '/storage/gjor/Data/kirchnerj/Data/curlyTab.csv';
        params.TRACESpath = '/storage/gjor/Data/kirchnerj/Data/Curly_TRACES/';
        %params.TRACESpath = '/storage/gjor/Data/kirchnerj/Data/Curly_TEST/';
        params.conditions = {'VISUAL','SENSORY','ROSTROLATERAL','ANTERIOLATERAL' , 'SUBPAR'};
        params.ExpNameStart = 12;
        params.plotRastl = 12;
        params.ID1s = 14;
        params.ID1e = 8;
        params.ID2s = 17;
        params.ID2e = 4;
        
        params.fixThresh = 95;
        params.diffEvRatePCT = 94;
        params.leaveOutRatio = 0.9;
        params.trDiffRat = 1.2;
        params.endEventPR = 1.7; % multiply event threshold by this number to get end threshold
        params.evRateWindow = 41;

        params.COND = 0;
        params.labels = {'V1' , 'S1' , 'RL' , 'AL'};

    elseif strcmp(DATASET , 'Palo')
        params.ROIpath = '/storage/gjor/Data/kirchnerj/outDirPalo/mat/';
        params.tabPath = '/storage/gjor/Data/kirchnerj/Data/Palo.csv';
        params.TRACESpath = '/storage/gjor/Data/kirchnerj/outDirPalo/';
        %params.TRACESpath = '/storage/gjor/Data/kirchnerj/Data/Friederike/TEST/';
        params.conditions = {'CONTROL','OXYTOCIN','CORTEXBUFFER', 'SUBPAR'};
        params.ExpNameStart = 1;
        params.plotRastl = 1;
        params.ID1s = 1;
        params.ID1e = 7;
        params.ID2s = 8;
        params.ID2e = 4;
        
        params.fixThresh = 93;
        params.diffEvRatePCT = 97;
        params.leaveOutRatio = 1.;
        params.trDiffRat = 1.2;
        params.endEventPR = 1.7; % multiply event threshold by this number to get end threshold
        params.evRateWindow = 41;

        params.COND = 0;
        params.labels = {'before OT' , 'after OT'};

    end
    % naming conventions
    params.outTableNames = {'rates' , 'durations' , 'durationsFWHM' ,...
        'startTimes' , 'amps' , 'realStartTimes' , 'jitter' , 'cond' ,...
        'F0','F0quant','age','id','id2','T'};
    params.pcaTableNames = {'pcaPerc','pcaPercAlt','age','id','id2','cond'};
    params.cpdTableNames = {'n1','n2','n3','n4'};
    
    % define neighborhoods for correlation analysis
    params.lNeigh = 5; % micrometer
    params.mNeigh = 50;
    params.umNeigh = 100; 
    params.uNeigh = 180;
    
    % how many principle components for 80 percent of the variance?
    params.PCAbound = 80;
    params.pcaSmoothWindowSize = 3;
        
    % define L and H events
    params.LOWERBOUND = 0.08;
    params.lEventLower = 0.2;
    params.lEventUpper = 0.8;
    % maximal number of ROIs
    params.upperN = 200;
    % if less then dudslim events, exclude ROI
    params.dudslim = 1;

    % paramter a for high/ lowpass filter
    params.FilterA = 0.2;%0.4; %0.2;
    % How many ROIs to exclude at the edges
    params.EdgePix = 20;
    params.idLim = 20;
    params.corrThresh = 80;

    params.nParRateWindow = 9; % naive participation rate window

    params.shuffleNum = 1000;   % number of repetitions for AR shuffle analysis
    
    params.ScatterPointSize = 15; % plotting scatter point size
end