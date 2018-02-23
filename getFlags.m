function flags = getFlags(DATASET)
    flags = struct();
    if strcmp(DATASET , 'Fried')
        flags.correctF0 = 1; % divide by F0
        flags.filterOut = 1; % remove certain percentage of bad ROIs
        flags.highFilt = 0; % apply high pass filter
        flags.lowFilt = 1; % apply low pass filter
        flags.changeAtCondChange = 1; % use different threshold after condition change
        flags.indivThresh = 0; % use individualized threshold per cell
        flags.ampF0Flag = 1; % use normalized traces
        flags.detrend = 1; % apply detrending
    elseif strcmp(DATASET , 'Curly')
        flags.correctF0 = 1;
        flags.filterOut = 1;
        flags.highFilt = 0;
        flags.lowFilt = 1;
        flags.changeAtCondChange = 1;
        flags.indivThresh = 0;
        flags.ampF0Flag = 1;
        flags.detrend = 1;
    elseif strcmp(DATASET , 'Palo')
        flags.correctF0 = 1;
        flags.filterOut = 1;
        flags.highFilt = 1;
        flags.lowFilt = 1;
        flags.changeAtCondChange = 0;
        flags.indivThresh = 1;
        flags.ampF0Flag = 0;
        flags.detrend = 1;
    end
    flags.onlyPrime = 0; % only use prime animals (indicated in experiment sheet)
    flags.skipNORMAL = 0; % skip NORMAL condition
    flags.skipBINOCULAR = 1; % ...
    flags.skipFORSKOLIN = 1;
    flags.skipVISUAL = 0;
    flags.skipSENSORY = 0;
    flags.skipANTERIOLATERAL = 0;
    flags.skipROSTROLATERAL = 0;
    flags.skipCONTROL= 0;
    flags.skipOXYTOCIN = 0;
    flags.skipCORTEXBUFFER = 1;
    
    flags.skipSUBPAR = 1; % exclude bad animals (failed motion correction etc) (indicated in experiment sheet)
    
    flags.plotVisible = 'on'; % display plots for the user
    flags.plotRast = 0; % display raster plot
    flags.corrFlag = 1; % compute correlation
    flags.cpdFlag = 0; % compute conditional probability distribution
    flags.scatterPos = 0; % scatter ROI positions
    flags.corrPlot = 0; % plot correlations
    flags.plotARMA = 1; % display AR plot
    flags.meanCorr = 0; % use mean value of correlations instead of cell individual
    flags.storeCorr = 0; % store correlations in file
    flags.addRegLine = 0; % add regression line to scatter plot
    
    flags.fixThreshPerExp = 1; % use same threshold for all experiments on one animal
    
    flags.removeDoubleCellROI = 1;
end