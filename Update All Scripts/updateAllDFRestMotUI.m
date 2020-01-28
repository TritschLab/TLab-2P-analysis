%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% mean DF/F during rest and during motion.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% Updated data structure, re-saved under its original name in the .mat file
% 
% Author: Jeffrey March, 2018

%% Main Code

[trials, pathname] = uigetfile('*.mat','MultiSelect','on');

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

cd(pathname)
totalFiles = 0;

roi = {};

for trial = 1:length(trials);
	load(trials{trial});                               
	totalFiles = totalFiles + 1
    
    roi{1} = data.dF1;
    roi{2} = data.dF2;
      
    for cellType = 1:2;
        
        roi{cellType}.DFRestMot = zeros(size(roi{cellType}.dF,1),2);
        DFCell = zeros(size(roi{cellType}.dF,1),2);
        
        for isMoving = [0, 1];
            roiData = data;
            
            if isMoving == 0;
                roiData.indOnsets = data.indOnsetsRest;
                roiData.indOffsets = data.indOffsetsRest;
            end

            totalTime = sum(roiData.indOffsets-roiData.indOnsets)/data.framerate;           
            
            for onset = 1:length(roiData.indOnsets)                    
                onsetMat = roi{cellType}.dF(:,roiData.indOnsets(onset):roiData.indOffsets(onset));                
                for nCell = 1:size(onsetMat,1)
                    DFCell(nCell, isMoving + 1) =  DFCell(nCell, isMoving + 1) + sum(onsetMat(nCell,:));
                end
            end    
            
            roi{cellType}.DFRestMot(:,isMoving + 1) = DFCell(:,isMoving + 1)/totalTime;
        end
        
    end
    
    data.dF1 = roi{1};
    data.dF2 = roi{2};
    
    save(trials{trial},'data')
end