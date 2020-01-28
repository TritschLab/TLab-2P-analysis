function [maskRegistered, centroids] = registeringROI(movNames, cellType, tform, movies)

% [maskRegistered, centroids] = registeringROI(movNames, cellType, tform, movies)
%
% Summary:  This function is run by roiReg and creates a mask of the 
% registered locations of all cells in a FOV. It functions by registering 
% the locations of all cells in a field of view using the same transform 
% that registers the corresponding tif file of the FOV to the tif file of 
% the reference FOV.
%
% Inputs:
% 
% 'movNames' - the names of all the tiff files
%
% 'cellType' - the cell types to be registered
%
% 'tform' - the transform required to register each FOV to the reference
%
% 'movies' - the tiff images of the FOVs
% 
% Outputs:
% 
% 'maskRegistered' - a binary map of locations of all the registered,
% active cells in a FOV.
% 
% 'centroids' - a matrix containing the centroid coordinates of each of the
% active cells.
%
% Author: Jeffrey March, 2018 


% cd('C:\MATLAB\Calcium Data\');
centroids = {};
maskRegistered = {};
for dataset = 1:length(movNames)
    datafile = [movNames{dataset}(1:end-3),'mat'];
    load(datafile);
            
        if cellType == 1
            temp = data.dF1;
        else if cellType == 2
                temp = data.dF2;
            else if cellType == 3
                    temp = data.dF3;
                end
            end
        end

        
        mask = zeros(512, 512);
        for i = 1:length(temp.roiList)
            tempROI = temp.roiList(i).indBody(~isnan(temp.roiList(i).indBody));
            mask(tempROI) = i;
%     coordinates = [ceil(tempROI/512), mod(tempROI-1,512)+1];
%     centroids(i,:) = mean(coordinates,1);
%     mask = insertText(mask,centroids(i,:),i);
        end

        if dataset > 1
            maskRegistered{dataset} = imwarp(mask,tform{dataset},'OutputView',imref2d(size(movies{1})));
        else
            maskRegistered{dataset} = mask;
        end
                       
        for i = 1:1:length(temp.roiList)
            tempROI = find(maskRegistered{dataset} == i);
            coordinates = [ceil(tempROI/512), mod(tempROI-1,512)+1];
            centroids{dataset}(i,:) = mean(coordinates,1);
        end
        
end

end
