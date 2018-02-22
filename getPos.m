function pos = getPos(csvPath , roiPath , plotFlag)
    if nargin < 3
        plotFlag = 0;
    end
    [pathTo , regionName , ~] = fileparts(csvPath);
    [ ~ , expName , ~] = fileparts(pathTo);
    fullRoiPath = fullfile(roiPath , expName , ['RoiSet ' regionName '.zip']);
    if ~exist(fullRoiPath,'file')
        pos = [];
        return
    end
    tab = ReadImageJROI(fullRoiPath);
    pos = zeros(length(tab),2);
    for ii = 1:length(tab)
        pos(ii,1) = tab{ii}.vnRectBounds(1);
        pos(ii,2) = tab{ii}.vnRectBounds(2);
    end
    if plotFlag
        figure(); suptitle('Positions of ROIS');
        scatter(pos(:,2),-pos(:,1))
    end
end