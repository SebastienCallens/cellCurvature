%%%% Quantify sheet sagging 
%%%% Sebastien Callens
clear; close all; clc

%% Input
substrate = 'Spheres';
fullName = 'example_data/20191204_DiffBD8S2Spheres_1_top.tif';

thrCoeff = 0.35;

numImgs = numel(imfinfo(fullName));
ImgInfo = imfinfo(fullName);
ppmu = ImgInfo.YResolution; % Pixels per micron
for i=1:numImgs
    Imgtemp = imread(fullName,i);
    Imgtemp = double(Imgtemp);
    Img(:,:,i) = Imgtemp;
end

%% Extract middle slice
mid = Img(:,:,round(numImgs/2));

%% Create z-projection
proj = sum(Img,3);
projOld = proj;
%% Normalize the images
mid = mid/mean(Img,'all');
proj = proj/mean(Img,'all');

%% Create masks and quantify cell anchoring
% Draw figures and create rois and masks
% Draw corner vertices all the way to the edges of the image
figure
colormap hot
imagesc(mid)
axis equal
axis tight
roi = drawpolygon('Color','w','LineWidth',1);
input('Draw mask. Press enter to continue...');
mask = uint16(createMask(roi));
% Convert to double
proj = double(proj);
mask = double(mask);
% Exclude the "dark" regions in the corners
thr = thrCoeff*mean(proj(:));
maskExclude = imbinarize(proj,thr);
maskFull = maskExclude.*mask;
% Apply mask to projections
proj_m = maskFull.*proj;
rho = mean(proj_m(maskFull>0),'all');

%% Binarize anchors
for i = 3:numImgs-8
    Imgtemp = Img(:,:,i);
    Imgtemp = maskFull.*Imgtemp;
    ImgtempBW = imbinarize(mat2gray(Imgtemp),0.1);
    ImgBW(:,:,i) = ImgtempBW;
end
anchorFraction = sum(ImgBW,'all')/(size(ImgBW,3)*sum(maskFull==1,'all')); 

%% Quantify amount of cell sagging using mask
[~,sheetIdx] = max(mask,[],1); % Find indixes of first occurences of 1 along every column of the mask
horIdx = mean([sheetIdx(1),sheetIdx(end)]); % reference line to measure sheet sagging (mean between left and right attachment point)
[lowest,lowIdx] = max(sheetIdx);
delta = lowest-horIdx;
delta_mu = delta/ppmu; 

%% Plot 
figure
colormap hot
imagesc(mid)
axis equal
axis tight
title('Middle top slices')
hold on
plot([0,size(mid,2)],[horIdx,horIdx],'w--','LineWidth',2)
plot([lowIdx,lowIdx],[horIdx,lowest],'w','LineWidth',2)
txtStr = strcat('\Delta=',num2str(round(delta_mu,2)),'\mum'); 
text(1.02*lowIdx,mean([horIdx,lowest]),txtStr,'Color','w');

%% display results:
disp(['Delta: ',num2str(delta_mu)]);
disp(['rho: ',num2str(rho)]);
disp(['Anchor fraction: ',num2str(anchorFraction)]);