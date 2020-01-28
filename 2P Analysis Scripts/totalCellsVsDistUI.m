%% Summary:
% 
% This script calculates and plots the fraction of SPNs vs. the total
% distance traveled for each acquisition across each condition.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'totalDistMat' - a cell array containing the total distance traveled in
% each acquisition for each condition.
% 
% 'totalFracCells' - a cell array containing the fraction of all SPNs of 
% each cell type that are analyzed for each acquisition of each condition
% 
% a plot of the fraction of all SPNs vs. the total distance traeled for
% each acquisition across each condition
% 
% Author: Jeffrey March, 2018

%% Actual Code

totalDistMat = cell(1,4);
totalFracCells = cell(1,4);

for conditions = 1:4
    conditionStrings = {'Pre','Acute','Chronic','L-dopa'};
        
    [trials, pathname] = selectFiles();
    
    totalFiles = 0;
    
    totalDistMat{conditions} = zeros(length(trials),1);
    totalFracCells{conditions} = zeros(length(trials),2);
    
    roi = cell(1,2);    
    
    for trial = 1:length(trials);
        load(trials{trial});
        totalFiles = totalFiles + 1
        
        roi{1} = data.dF1;
        roi{2} = data.dF2;
        
        totalDistMat{conditions}(trial,1) = data.totalDistance;
        for cellType = 1:2
            totalFracCells{conditions}(trial,cellType) = size(roi{cellType}.dF,1)/roi{cellType}.totalSPNs;           
        end        
    end    

end


%%

for conditions = 1:4
    for cellType = 1:2
        figure(cellType);
        plot(totalDistMat{conditions},totalFracCells{conditions}(:,cellType)','x');
        hold on;
    end
end

