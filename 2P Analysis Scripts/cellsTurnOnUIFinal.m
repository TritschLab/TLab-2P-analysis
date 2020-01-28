%% Summary:
% 
% This script calculates and plots the cumulative distribution of the times 
% when each cell is first activated relative to distance traveled in a 
% given trial.
% 
% Properties for paper: 
% 
% conditions = 1 (this is for getting individual cumulative distributions
% for each FOV in each separate experimental condition)
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'finalMat' - cell array containing a matrix for each cell type, with the 
% first row being the distance where each new cell turns on and the second 
% row being the number cell that is in that FOV
% 
% 'totalDistances' - cell array containing a matrix for each experimental
% condition with the total distance traveled for each FOV
% 
% 'indivCumDistrib' - cell array containing a matrix for each cell type
% with the individual cumulative distribution for each FOV, the last
% non-NaN value of each row being the total distance traveled, NOT the last
% cell turning on.
% 
% plot of individual cumulative distributions for each FOV
% 
% Author: Jeffrey March, 2018

%%
totalDistances = cell(1,4);

% Iterating through each condition
for conditions = 1 %1:3
    conditionStrings = {'Pre','Acute','Chronic','L-dopa'};
        
    [trials, pathname] = uigetfile('*.mat','MultiSelect','on');
    
    % In case you want to look at one acquisition
    if ~iscell(trials)
        tempTrials = trials;
        trials = cell(1);
        trials{1} = tempTrials;
    end
    
    cd(pathname)
    
    % Initializing values
    totalFiles = 0;
    npeaks = 1; % Only counting the first spike of each cell
    MINW = 4;
    noiseThresh = 0.4;
    stdAway = 5;
    roi = cell(1,2);
    warning('off');
    finalMat = cell(1,2);
    
       
    for trial = 1:length(trials);
        load(trials{trial});
        totalFiles = totalFiles + 1
        
        roi{1} = data.dF1;
        roi{2} = data.dF2;
        
        for cellType = 1:2;
            LOCS = zeros(1,size(roi{cellType}.dF,1));
            PKS = zeros(1,size(roi{cellType}.dF,1));     
            
            roi2 = zeros(size(roi{cellType}.dF));
            % Smoothing Cells
            for cellNum = 1:size(roi{cellType}.dF,1)
                roi2(cellNum,:) = smooth(roi{cellType}.dF(cellNum,:),5)';
            end
            
            % Finding all peaks in all cells
            for cellNum = 1:size(roi2,1)
                MPP = mean(roi2(cellNum,roi2(cellNum,:) < noiseThresh)) + stdAway*std(roi2(cellNum,roi2(cellNum,:) < noiseThresh));
                MPH = MPP;
                if ~isempty(findpeaks(roi2(cellNum,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW,'NPeaks',npeaks))
                    [PKS(cellNum), LOCS(cellNum)] = findpeaks(roi2(cellNum,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW,'NPeaks',npeaks);
                end
            end
            
            LOCS = LOCS(find(LOCS ~= 0));
            
            distForCell = zeros(2,length(LOCS));
            distForCell(2,:) = (1:length(LOCS)); %/length(LOCS); %add /length(LOCS) to make y axis "fraction of active cells"
            cumDist = cumsum(abs(data.vel))/data.framerate;
            distForCell(1,:) = sort(cumDist(LOCS)); %/cumDist(end); %add /cumDist(end) to make x axis "fraction of total distance"
            
            finalMat{cellType} = cat(2,finalMat{cellType},distForCell);
        end
        
        if data.totalDistance > 50
            display([trials{trial},' ',num2str(data.totalDistance)]);
        end
        
        totalDistances{conditions} = [totalDistances{conditions}, data.totalDistance];
        
    end
    
%     %Run this for all conditions, one after the other, to get full picture
%     for cellType = 1:2
%         figure(cellType + 40);
%         plot(finalMat{cellType}(1,:),finalMat{cellType}(2,:),'x')
%         title(['Fraction of Active D',num2str(cellType),' Cells Activated vs. Distance Traveled']);
%         xlabel('Distance Traveled (m)');
%         ylabel('Fraction of Active Cells Activated');
%         legend(conditionStrings{1},conditionStrings{2},conditionStrings{3},conditionStrings{4});
%         %     xlim([0,1]);
% %         ylim([0,1]);
%         hold on
%         
%     end
    
%     for cellType = 1:2
%         [dist, distSortOrder] = sort(finalMat{cellType}(1,:));
%         cellNum = finalMat{cellType}(2,distSortOrder);
%         
%         k = lsqcurvefit(@logistic,3,dist,cellNum)
%         
%         fitLine = logistic(k, dist);
%         
%         figure(cellType + 10);
%         hold on
%         plot(dist,fitLine);
%         title(['Fraction of Active D',num2str(cellType),' Cells Activated vs. Distance Traveled (Fit)']);
%         xlabel('Distance Traveled (m)');
%         ylabel('Fraction of Active Cells Activated');
% %         ylim([0,1]);
%         legend(conditionStrings{1},conditionStrings{2},conditionStrings{3},conditionStrings{4});
%     end

    % Calculating distance between cells turning on
%     distBetweenCells = cell(1,2);
%     for cellType = 1:2
%         distBetweenCells{cellType} = nan(sum(finalMat{cellType}(2,:) == 1),max(finalMat{cellType}(2,:)));
%         row = 0;
%         for i = 1:size(finalMat{cellType},2)
%             if finalMat{cellType}(2,i) == 1
%                 row = row + 1;
%                 distBetweenCells{cellType}(row,finalMat{cellType}(2,i)) = finalMat{cellType}(1,i);
%             else
%                 distBetweenCells{cellType}(row,finalMat{cellType}(2,i)) = finalMat{cellType}(1,i) - finalMat{cellType}(1,i-1);
%             end
%         end
% 
%         figure;
%         plot(smooth(nanmean(distBetweenCells{cellType},1),3)');
%         xlim([0,120]);
%         ylim([0,20]);
%     end
    
end


%% Creating individual cumulative distributions

indivCumDistrib = cell(1,2);
for cellType = 1:2
    indivCumDistrib{cellType} = nan(sum(finalMat{cellType}(2,:) == 1),1+max(finalMat{cellType}(2,:)));
    row = 0;
    for i = 1:size(finalMat{cellType},2)
        if finalMat{cellType}(2,i) == 1
            if i > 1
                indivCumDistrib{cellType}(row,finalMat{cellType}(2,i-1) + 1) = totalDistances{1}(row);
            end
            row = row + 1;
        end
        indivCumDistrib{cellType}(row,finalMat{cellType}(2,i)) = finalMat{cellType}(1,i);
    end
    indivCumDistrib{cellType}(row,finalMat{cellType}(2,i) + 1) = totalDistances{1}(row);
    
    
    for row = 1:size(indivCumDistrib{cellType},1)
        if sum(isnan(indivCumDistrib{cellType}(row,:))) > 0
            numCells = find(isnan(indivCumDistrib{cellType}(row,:)),1) - 2;
        else
            numCells = length(indivCumDistrib{cellType}(row,:)) - 1;
        end
        if numCells > 0
            percentOfTotal = (1/numCells):(1/numCells):1;
            figure(50 + cellType);
            hold on
            prev = plot([0,indivCumDistrib{cellType}(row,1:numCells)],[0,percentOfTotal],'-');
            color = get(prev, 'Color');
            plot(indivCumDistrib{cellType}(row,numCells:numCells+1),[1,1 + 1/numCells],'--','Color',color);
            hold on
        end
    end
end


