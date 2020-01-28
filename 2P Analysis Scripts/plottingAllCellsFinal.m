%% Summary:
% 
% This script plots all cells in a FOV in various ways: plotting all cell
% traces vertically, with the wheel beneath and movement bouts hightlighted; 
% as a heat plot where each row is a cell; as a raster plot with movement 
% bouts highlighted; and as an average DF/F trace with movement bouts 
% hightlighted.
% 
% Properties for paper: 
% 
% SW024_170502_Acq3 (mouse used for paper)
% 
% Inputs:
% 
% This will operate on the open 'data' file
%
% Outputs:
% 
% plot of cells, dSPNs in red, iSPNs in green, layered vertically above the
% wheel trace, with movement being highlighted
% 
% smoothed heatmap of cells, dSPNs and iSPNs in different heatmaps
% 
% raster plot of cells, dSPN spikes during movement in red, iSPN spikes 
% during movement in green, and all SPN spikes during rest in black, with 
% movement being highlighted 
% 
% average DF/F traces, with dSPN average in red and iSPN average in green,
% with movement bouts highlighted
% 
% Author: Jeffrey March, 2018


%% Plotting all Cells Calcium w/Wheel + imagesc
dist = 4;
offset = -5;
scale = 20;
figure;
hold on;
totalCells = cat(1,scale*data.vel,data.dF1.dF,data.dF2.dF); %data.dF2.dF(1:end - 4,:) is used for the paper this acquisition

numD1 = size(data.dF1.dF,1);
threshold = .002;
mov = abs(data.vel) > threshold;
for i = 1:size(totalCells,1)
   
    if i == 1
        plot(data.frameTime,scale*data.vel - 8);
    else if i <= numD1 + 1
            plot(data.frameTime,offset + totalCells(i,:),'Color',[0.8,0,0]);
%         plot(data.frameTime,offset + atRest,'r.','MarkerSize',0.5);
%         plot(data.frameTime,offset + inMot,'k.','MarkerSize',0.5);
        else
            plot(data.frameTime,offset + totalCells(i,:),'Color',[0,0.8,0]);
%         plot(data.frameTime,offset + atRest,'g.','MarkerSize',0.5);
%         plot(data.frameTime,offset + inMot,'k.','MarkerSize',0.5);
        end
    end
    offset = offset + dist;
end
a = area(data.frameTime,dist*size(totalCells,1)*mov,'FaceColor','k','EdgeColor','none'); % Shading movement bouts
a.FaceAlpha = 0.1;

% Smoothing matrix for heat map
smoothedCells = zeros(size(totalCells));
smoothing = round(2*data.framerate);
for i = 1:size(totalCells,1)
%     smoothedCells(i,:) = totalCells(i,:)/max(totalCells(i,:)); % Normalizing signal
    smoothedCells(i,:) = smooth(totalCells(i,:),smoothing)';
end
% imagesc(smoothedCells,[0,.5]);

figure;
colormap('bone');
colorLimits = [.15, .3];
subplot(3,1,1);
imagesc(totalCells(1,:)) % Plotting wheel

subplot(3,1,2);
imagesc(smoothedCells(2:numD1+1,:),colorLimits); % Plotting dSPNs

subplot(3,1,3);
imagesc(smoothedCells(numD1+1:end,:),colorLimits); % Plotting iSPNs

%% Plotting Raster

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi{1} = data.dF1;
roi{2} = data.dF2;
counter = 0;
figure;
title('Calcium Spikes, Raster Plot')
ylabel('Cell #')
xlabel('Time (s)')
hold on
threshold = .002;
shift = 0;
mov = abs(data.vel) > threshold;
for iShift = 1:shift
    mov = [mov(1+shift:end),zeros(1,shift)] | mov | [zeros(1,shift),mov(1:end-shift)];
end

for cellType = 1:2;
    
    roi{cellType}.peakFreqRestMot = zeros(size(roi{cellType}.dF,1),2);
    numSpikesCell = zeros(size(roi{cellType}.dF,1),2);
    
    roiData = data;
    
    roi2 = zeros(size(roi{cellType}.dF));
    % Z-scoring Cells
    for nCell = 1:size(roi{cellType}.dF,1)
        roi2(nCell,:) = smooth(roi{cellType}.dF(nCell,:),5)';
    end
    
    for cellROI = 1:size(roi2,1)
        inMot = NaN(1,size(roi2,2));
        atRest = NaN(1,size(roi2,2));
        raster = NaN(1,size(roi2,2));
        warning('OFF');
        
        MPP = mean(roi2(cellROI,roi2(cellROI,:) < noiseThresh)) + stdAway*std(roi2(cellROI,roi2(cellROI,:) < noiseThresh));
        baseline = mean(roi2(cellROI,roi2(cellROI,:) < noiseThresh));
        MPH = MPP;
        [tempPKS, tempLOCS] = findpeaks(roi2(cellROI,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
        inMot(tempLOCS(find(mov(tempLOCS)))) = counter;
        atRest(tempLOCS(find(~mov(tempLOCS)))) = counter;
        raster(tempLOCS) = counter;
        if cellType == 1
            plot(data.frameTime,atRest,'k.') % Black when spike at rest
            plot(data.frameTime,inMot,'r.') % Red when spike at movement
        else
            plot(data.frameTime,atRest,'k.') % Black when spike at rest
            plot(data.frameTime,inMot,'g.') % Green when spike at movement
        end
 
        counter = counter + 1;
    end

end

a = area(data.frameTime,counter*mov,'FaceColor','k','EdgeColor','none');
a.FaceAlpha = 0.1;

%% Plotting Vel and Mean DF

figure;
subplot(2,1,1);
plot(data.frameTime,smooth(data.vel,2*data.framerate)');
subplot(2,1,2);
hold on
plot(data.frameTime,smooth(mean(data.dF1.dF,1),2*data.framerate)','r');
plot(data.frameTime,smooth(mean(data.dF2.dF,1),2*data.framerate)','g');
threshold = .002;
mov = abs(data.vel) > threshold;
% plot(downsample(data.frameTime,5), downsample(.16*mov,5),'k','LineWidth',.1);
a = area(downsample(data.frameTime,5),downsample(.16*mov,5),'FaceColor','k','EdgeColor','none');
a.FaceAlpha = 0.1;
