%%%% Quantify RUNX2 intensity in a series of ROIs for a set of images
%%%% Sebastien Callens

clear; close all; clc
cmap = coolwarm(256);
color1 = [239,193,100]/255; % yellowish
color2 = [243,131,93]/255; % orangeish
color3 = [243,89,85]/255; % redish
color4 = [40,98,117]/255; % dark blueish
color5 = [0,67,76]/255; % darker blueish
%% Load images
% Select substrate and experiment case and load set of image names
substrate = 'Unduloid';
caseExp = 'D8ConvexDiff';
imgName = selectImgName(substrate,caseExp);
% load curvature maps
load(strcat(substrate,'_k1Map.mat'));
load(strcat(substrate,'_k2Map.mat'));
gauss = curvMap_k1.*curvMap_k2;

%% select ROI centres
switch substrate
    case 'Unduloid'
        xval = 297;
        yval = [505,928,1351,1774,2197,2620,3043,3466];
end
ROIsize = 110; % Size of ROI in pixels

%% Loop over all images and ROIs
for i = 1:length(imgName)
    runx2_name = strcat(imgName{i}(1:end-5),'runx2');
    dapi_name = strcat(imgName{i}(1:end-5),'nucleus');
    runx2 = load(strcat(runx2_name,'.mat'));
    runx2 = runx2.Img3;
    dapi = load(strcat(dapi_name,'.mat'));
    dapi = dapi.Img2;
    runx2_store(:,:,i) = runx2;
    dapi_store(:,:,i) = dapi;
    % Run over all ROIs
    for j = 1:length(yval)
        xcoord = [xval-ROIsize/2,xval+ROIsize/2,xval+ROIsize/2,xval-ROIsize/2];
        ycoord = [yval(j)-ROIsize/2,yval(j)-ROIsize/2,yval(j)+ROIsize/2,yval(j)+ROIsize/2];
        % Mask images
        cropMask = roipoly(runx2,xcoord,ycoord);
        % Now remove black pixels that result after applying the mask
        runx2ROI = cropMasked(runx2,cropMask);
        dapiROI = cropMasked(dapi,cropMask);
        gaussROI = cropMasked(gauss,cropMask);
        % Gauss filter dapi ROI
        dapiROI_gauss = imgaussfilt(dapiROI,2);
        % binarize
        dapiROI_bin = imbinarize(dapiROI_gauss);
        % Apply mask to runx2 and dapi ROIs
        % Apply binarized image as mask
        runx2_masked = runx2ROI.*uint8(dapiROI_bin);
        dapi_masked = dapiROI.*uint8(dapiROI_bin);
        gauss_masked = gaussROI.*double(dapiROI_bin);
        % calculate mean intensity of runx2 and dapi
        runx2_int((i-1)*length(yval)+j) = mean(runx2_masked(runx2_masked>0));
        dapi_int((i-1)*length(yval)+j) = mean(dapi_masked(dapi_masked>0));
        gauss_mean((i-1)*length(yval)+j) = mean(gauss_masked,'all');
        % store runx2ROI, dapiROI, dapiROI_bin
        runx2ROI_store((i-1)*length(yval)+j,:,:) = runx2ROI;
        dapiROI_store((i-1)*length(yval)+j,:,:) = dapiROI;
        dapiROI_bin_store((i-1)*length(yval)+j,:,:) = dapiROI_bin;
    end
    % Normalize the last set of ROIs by the mean of that set
    runx2_int((i-1)*length(yval)+1:i*length(yval)) = runx2_int((i-1)*length(yval)+1:i*length(yval))/mean(runx2_int((i-1)*length(yval)+1:i*length(yval)));
    dapi_int((i-1)*length(yval)+1:i*length(yval)) = dapi_int((i-1)*length(yval)+1:i*length(yval))/mean(dapi_int((i-1)*length(yval)+1:i*length(yval)));
end
        
%% Store data
% create vectors with runx2 intensities at the pos and neg Gaussian
% curvatures
runx2_pos = runx2_int(gauss_mean>0);
runx2_neg = runx2_int(gauss_mean<0);
dapi_pos = dapi_int(gauss_mean>0);
dapi_neg = dapi_int(gauss_mean<0);
% Fill out with NaNs so that they can be pasted in the data matrix
runx2_pos(end+1:length(runx2_int)) = NaN;
runx2_neg(end+1:length(runx2_int)) = NaN;
dapi_pos(end+1:length(runx2_int)) = NaN;
dapi_neg(end+1:length(runx2_int)) = NaN;
% Create data matrix
% intData = [runx2_int',dapi_int',gauss_mean',runx2_pos',runx2_neg',dapi_pos',dapi_neg'];
% saveStr = strcat('Data/RUNX2_',caseExp,'_',substrate,'.xlsx');
% xlswrite(saveStr,{'RUNX2','DAPI','Gauss','RUNX2_K>0','RUNX2_K<0','DAPI_K>0','DAPI_K<0'},'Sheet1','A1');
% xlswrite(saveStr,intData,'Sheet1','A2');

%% Verify ROIs on example image
figure;
imagesc(runx2); axis equal; axis tight; colormap(cmap); colorbar; set(gca,'visible','off')
hold on
for j = 1:length(yval)
    xcoord = [xval-ROIsize/2,xval+ROIsize/2,xval+ROIsize/2,xval-ROIsize/2];
    ycoord = [yval(j)-ROIsize/2,yval(j)-ROIsize/2,yval(j)+ROIsize/2,yval(j)+ROIsize/2];
    patch('XData',xcoord,'YData',ycoord,'edgecolor','k','facecolor','none','Linewidth',1);
end