%% Summary:
% 
% This script performs calculations on cells that have previously been
% registered and tracked across experimental days. Nothing came of this in
% the paper, but there is the backbone for analysis of all peak stats
% across days. To start, follow the instructions in 'Setup' below. This
% code should be refined and adapted.
% 
% Inputs:
% 
% user must define allAcqs from excel template
%
% Outputs:
% 
% 'freqRest10' - frequency during rest for each cell
% 
% 'freqMov10' - frequency during movement for each cell
% 
% 'height10', 'rise10', 'decay10' - peak stats for each cell
% 
% 'identity' - the cell type (determined probabilistically) of each cell
% 
% 'restMovActive10' - is the cell active during rest, motion, both, or neither
% 
% 'motorSpontDiffRest', 'motorSpontDiffMov' - difference between motor and 
% spontaneous spike frequencies during rest and movement
% 
% Author: Jeffrey March, 2018


%% Setup

% Instructions:
allAcqs = {};
% Copy and paste template from excel sheet into allAcqs. Run first three
% sections of this code.


%% Creating Freq, Height, Rise, and Decay Stats Matrix

cd('L:\tritsn01labspace\Marta\2p_Data\Calcium Data\')
freqRest = zeros(size(allAcqs,1)-1,size(allAcqs,2));
freqMov = zeros(size(allAcqs,1)-1,size(allAcqs,2));
height = zeros(size(allAcqs,1)-1,size(allAcqs,2));
rise = zeros(size(allAcqs,1)-1,size(allAcqs,2));
decay = zeros(size(allAcqs,1)-1,size(allAcqs,2));
identity = zeros(size(allAcqs,1)-1,1);
motor = 1:4;
spont = 6:9;

for j = 1:size(allAcqs,2)
    if exist(allAcqs{1,j})
        load(allAcqs{1,j})
        roi = {data.dF1,data.dF2, data.dF3};
    end   
    for i = 2:size(allAcqs,1)
        if allAcqs{i,j} > 0
            if floor(allAcqs{i,j}/1000) < 3
                freqRest(i-1,j) = roi{floor(allAcqs{i,j}/1000)}.peakFreqRestMot(mod(allAcqs{i,j},1000),1);
                freqMov(i-1,j) = roi{floor(allAcqs{i,j}/1000)}.peakFreqRestMot(mod(allAcqs{i,j},1000),2);
                height(i-1,j) = roi{floor(allAcqs{i,j}/1000)}.peakHeights(mod(allAcqs{i,j},1000));
                rise(i-1,j) = roi{floor(allAcqs{i,j}/1000)}.riseTimes(mod(allAcqs{i,j},1000));
                decay(i-1,j) = roi{floor(allAcqs{i,j}/1000)}.decayTimes(mod(allAcqs{i,j},1000));
                identity(i-1) = identity(i-1) + 2*(floor(allAcqs{i,j}/1000) - 1.5);
            else
                identity(i-1) = nan;
            end
        else
            if ~isempty(allAcqs{i,j})
                freqRest(i-1,j) = allAcqs{i,j};
                freqMov(i-1,j) = allAcqs{i,j};
            else
                freqRest(i-1,j) = nan;
                freqMov(i-1,j) = nan;
            end
            height(i-1,j) = nan;
            rise(i-1,j) = nan;
            decay(i-1,j) = nan;
        end

    end
end

identity(find(identity > 0)) = 2;
identity(find(identity < 0)) = 1;
identity(find(identity == 0)) = nan;

%% Creating matrices of 10 conditions, separating D1 and D2, plotting

freqRest10 = {zeros(length(find(identity == 1)),10), zeros(length(find(identity == 2)),10)};
freqMov10 = {zeros(length(find(identity == 1)),10), zeros(length(find(identity == 2)),10)};
height10 = {zeros(length(find(identity == 1)),10), zeros(length(find(identity == 2)),10)};
rise10 = {zeros(length(find(identity == 1)),10), zeros(length(find(identity == 2)),10)};
decay10 = {zeros(length(find(identity == 1)),10), zeros(length(find(identity == 2)),10)};

indCellType = [1,1];
for i = 1:size(freqMov,1)
    if identity(i) == 1 || identity(i) == 2
        for j = 1:10
            freqRest10{identity(i)}(indCellType(identity(i)),j) = freqRest(i,j*3 - 2);
            freqMov10{identity(i)}(indCellType(identity(i)),j) = freqMov(i,j*3 - 2);
            height10{identity(i)}(indCellType(identity(i)),j) = height(i,j*3 - 2);
            rise10{identity(i)}(indCellType(identity(i)),j) = rise(i,j*3 - 2);
            decay10{identity(i)}(indCellType(identity(i)),j) = decay(i,j*3 - 2);           
        end
        indCellType(identity(i)) = indCellType(identity(i)) + 1;
    end
end

% for cellType = 1:2
%     figure;
%     hold on
%     condits = spont;
%     for i = 1:size(freqRest10{cellType},1)
%         plot(freqRest10{cellType}(i,condits),'x-')
% %         plot(freqMov10{cellType}(i,condits),'o-')
%     end
% end

clearvars i j roi
%% Classifying cells as active during rest or during movement

restMovActive10 = {};
for cellType = 1:2
    restMovActive10{cellType} = NaN(size(freqMov10{cellType})); % default is NaN
    restMovActive10{cellType}(find(freqMov10{cellType} == 0 & freqRest10{cellType} == 0)) = 1; % active during neither
    restMovActive10{cellType}(find(freqMov10{cellType} == 0 & freqRest10{cellType} > 0)) = 2; % active during rest only
    restMovActive10{cellType}(find(freqMov10{cellType} > 0 & freqRest10{cellType} == 0)) = 3; % active during movement only
    restMovActive10{cellType}(find(freqMov10{cellType} > 0 & freqRest10{cellType} > 0)) = 4; % active during both
end

%% Comparing heights across trials

%% Comparing motor vs. spontaneous stats

numCondits = 5;
motorSpontDiffRest = {};
motorSpontDiffMov = {};
for cellType = 1:2
    resultRest = freqRest10{cellType}(:,1:numCondits) - freqRest10{cellType}(:,1+numCondits:end);
    motorSpontDiffRest{cellType} = resultRest(find(~isnan(resultRest)));
    resultMov = freqMov10{cellType}(:,1:numCondits) - freqMov10{cellType}(:,1+numCondits:end);
    motorSpontDiffMov{cellType} = resultMov(find(~isnan(resultMov)));

    figure;
    histogram(resultMov,100);
    title(['D',num2str(cellType),' freqMov Motor vs. Spont STDev: ',num2str(std(motorSpontDiffMov{cellType}))]);
    figure;
    histogram(resultRest,100);
    title(['D',num2str(cellType),' freqRest Motor vs. Spont STDev: ',num2str(std(motorSpontDiffRest{cellType}))]);
end
%% Comparing pre and post lesion motor

condits = spont;
for cellType = 1:2
    numCells = 0;
    figure;
    meanSignalRest = zeros(1,length(condits));
    meanSignalMov = zeros(1,length(condits));
    for i = 1:size(freqMov10{cellType},1)    
        if ~isnan(mean(freqMov10{cellType}(i,condits)))
            meanSignalRest = meanSignalRest + freqRest10{cellType}(i,condits);
            meanSignalMov = meanSignalMov + freqMov10{cellType}(i,condits);
            numCells = numCells + 1;
            
            ax1 = subplot(2,1,1);
            plot(freqRest10{cellType}(i,condits),'x-')
            title(['D',num2str(cellType),' Rest Freq Across Conditions'])
            hold on
            
            ax2 = subplot(2,1,2);
            plot(freqMov10{cellType}(i,condits),'x-')
            title(['D',num2str(cellType),' Movement Freq Across Conditions'])
            hold on
        end       
    end
    
    ax1 = subplot(2,1,1);
    plot(meanSignalRest/numCells,'r^-')
    ylim([0,0.11]);
    ax2 = subplot(2,1,2);
    plot(meanSignalMov/numCells,'r^-')
    ylim([0,0.11]);
    
    linkaxes([ax1,ax2],'x')
end
%% Calculating number of cells that are new/the same

% Cell has 2 columns (cell types) and 4 rows (cells appear, cells remain,
% cells disappear, total cells).
 
condits = spont;
cellsAddedSame = cell(4,2);

for cellType = 1:2
    for col = condits(1:end-1)
        cellsAddedSame{1,cellType}(col - condits(1) + 1) = sum(isnan(freqMov10{cellType}(:,col)) & ~isnan(freqMov10{cellType}(:,col + 1)));
        cellsAddedSame{2,cellType}(col - condits(1) + 1) = sum(~isnan(freqMov10{cellType}(:,col)) & ~isnan(freqMov10{cellType}(:,col + 1)));
        cellsAddedSame{3,cellType}(col - condits(1) + 1) = sum(~isnan(freqMov10{cellType}(:,col)) & isnan(freqMov10{cellType}(:,col + 1)));
        cellsAddedSame{4,cellType}(col - condits(1) + 1) = sum(~isnan(freqMov10{cellType}(:,col)));
    end
end

transitions = zeros(4,3);
for cellType = 1:2
    for i = 1:3
        transitions(cellType*2 - 1,i) = cellsAddedSame{1,cellType}(i)/cellsAddedSame{4,cellType}(i);
        transitions(cellType*2,i) = cellsAddedSame{3,cellType}(i)/cellsAddedSame{4,cellType}(i);
    end
end


%% Calculating number of cells that are different from prelesion

% Cell has 2 columns (cell types) and 4 rows (cells appear, cells remain,
% cells disappear, total cells).

condits = spont;
cellsAddedSame = cell(4,2);

for cellType = 1:2
    for col = condits(2:end)
        cellsAddedSame{1,cellType}(col - condits(1)) = sum(isnan(freqMov10{cellType}(:,1)) & ~isnan(freqMov10{cellType}(:,col)));
        cellsAddedSame{2,cellType}(col - condits(1)) = sum(~isnan(freqMov10{cellType}(:,1)) & ~isnan(freqMov10{cellType}(:,col)));
        cellsAddedSame{3,cellType}(col - condits(1)) = sum(~isnan(freqMov10{cellType}(:,col)));
        cellsAddedSame{4,cellType}(col - condits(1)) = sum(~isnan(freqMov10{cellType}(:,1)));
    end
end

transitions = zeros(4,3);
for cellType = 1:2
    for i = 1:3
        transitions(cellType*2 - 1,i) = cellsAddedSame{1,cellType}(i)/cellsAddedSame{4,cellType}(i);
        transitions(cellType*2,i) = cellsAddedSame{2,cellType}(i)/cellsAddedSame{4,cellType}(i);
    end
end
%%
% 
% for i = 1:size(active,1)
%     for j = 1:size(active,2)
%         %finding total active acquisitions for each cell
%         if ~isempty(height{i,j})
%             activeAcqs(i) = activeAcqs(i) + 1;
%             activeCondit(ceil(j/3)) = activeCondit(ceil(j/3)) + 1;
%         end
%         
%         %finding total active pairs
%         if mod(j-1,6) == 0 && active(i,j) && active(i,j+3)
%             activePairs(i) = activePairs(i) + 1;
%         end
%         
%         %finding total active motor acquisitions for each cell
%         if active(i,j) && j <= 15
%             activeMot(i) = activeMot(i) + 1;
%         end
%         
%         %finding total active spontaneous acquisitions for each cell
%         if active(i,j) && j > 15
%             activeSpont(i) = activeSpont(i) + 1;
%         end
%     end
% end
% 
% figure;
% hold on
% % plot(numActive,'b')
% % plot(activePairs, 'r')
% plot(activeMot, 'g+')
% plot(activeSpont, 'ko')
% 
% figure;
% hold on
% plot(activeCondit(1:5),'b')
% plot(activeCondit(6:10),'r')
