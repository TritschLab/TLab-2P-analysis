function MotC_AVGs(folder,Motif,maxShift,ref)
%now we apply all the average motion correction to the averages
% then we go and put that on each motif

disp('Motion Correcting the Averages')
%motion correct the averaged tiff
avg_name=[folder,'4-Avgs\AVG.tif'];
is = imageSeries(avg_name);
[results,~] = is.motionCorrect('savePath',[folder,'4-Avgs\AVGMotC.tif'],'referenceFrame',ref,'maxShift',maxShift);
xshifts=results.xshifts;
yshifts=results.yshifts;
save([folder,'4-Avgs\MakeTiffparams.mat'],'xshifts','yshifts')


newT_fold=[folder,'5-FinalMotifs\'];
if ~exist(newT_fold,'dir')
    mkdir(newT_fold)
end


%based off of how we are correcting, what portion should we keep?
%maximum cut offs
maxX=ceil(max(xshifts));%cut off on the left
maxY=ceil(max(yshifts));%cut off on top
minX=floor(min(xshifts));%cut off on the right
minY=floor(min(yshifts));%cut off the bottom

%current size
origT_name=[folder,'4-Avgs\AVG.tif'];%the name/location of the tiff file
origT_obj = Tiff(origT_name,'r');%create a tiff object, only for reading
l=origT_obj.getTag('ImageLength');%rows (first index)
w=origT_obj.getTag('ImageWidth');%columns (second index). i.e. it would be [len,wid]=sizeimage)
origT_obj.close();

%what to keep
keepRangeY = max(1,1+maxY):min(l,l+minY);
keepRangeX = max(1,1+maxX):min(w,w+minX);
keepRangeY=keepRangeY(1:end-1);% for some reason they keep on cutting off a part of it

motifNames={Motif(:).name};
indReport=floor((.1:.1:1)*(length(motifNames)));
fprintf('Applying to Motifs: ')
for m=1:length(motifNames)
%     clear Allavg1
    origT_name=[folder,filesep,'3a-MotifsMin\',filesep,motifNames{m}];%the name/location of the tiff file
    origT_obj = Tiff(origT_name,'r');%create a tiff object, only for reading
    x=xshifts(m);%this is this files shift relative to the zeroth file
    y=yshifts(m);
    newT_name=[newT_fold,motifNames{m}];
    numF=diff(Motif(m).frames)+1;
    for i=1:numF
        origT_obj.setDirectory(i);
        imageData = origT_obj.read();%get the image
        %remember, x and y are flipped. y is row, x is column
        X=(1:l)+y;        Y=((1:w)+x)';
        Xq=1:l;        Yq=(1:w)';
        %remember they flip x and y
        imageData2=interp2(Y,X,double(imageData),Yq,Xq,'linear');%interpolate
        imageData2=uint16(imageData2(keepRangeY,keepRangeX));%cut out the edges
        if i==1 
            imwrite(imageData2,newT_name,'TIF','compression','none')
        else
            imwrite(imageData2,newT_name,'TIF','WriteMode','append','compression','none')
        end
    end
    if sum(indReport==m)
        fprintf([num2str(m/length(motifNames),1),','])
    end
end
disp('Done')