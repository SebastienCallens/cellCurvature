function [Data_1D] = intensityVSCurvatureFun(substrate,imgName,caseExp,pixSize)
% This function is used to batch-process sets of images and store the 1D
% data, to be processed in Prism later on
%% Load processed images
for i = 1:length(imgName)
    Imgtemp = load(strcat(imgName{i},'.mat'));
    Imgtemp = Imgtemp.Img;
    % Convert image to double
    Imgtemp = double(Imgtemp);
    Imgtemp_norm = Imgtemp/mean(Imgtemp,'all');
    % Store non-normalized and normalized versions
    ImgIn(:,:,i) = Imgtemp;
    ImgIn_norm(:,:,i) = Imgtemp_norm;
end
Img = ImgIn_norm;
%% Load curvature map
AA = load(strcat(substrate,'_k1Map.mat'));
BB = load(strcat(substrate,'_k2Map.mat'));
curvMap_k1 = AA.curvMap_k1;
curvMap_k2 = BB.curvMap_k2;
if strcmp(substrate,'Wavy')
    curvMap_k1 = flip(curvMap_k1,1);
    curvMap_k2 = flip(curvMap_k2,1);
end
% Normalize curvature by radius of sphere
curvMap_k1 = curvMap_k1*180;
curvMap_k2 = curvMap_k2*180;
if strcmp(caseExp(3:9),'Concave')
    k1temp = -1*curvMap_k1;
    k2temp = -1*curvMap_k2;
    curvMap_k1 = k2temp;
    curvMap_k2 = k1temp;
end

gaussCurv = curvMap_k1.*curvMap_k2;
meanCurv = 0.5*(curvMap_k1+curvMap_k2);
%% Create Euclidean distance transform for the k2 curvature maps
k2_binA = imbinarize(curvMap_k2,-0.00000001);
k2_binB = k2_binA;
k2_binA(curvMap_k1==0) = 0;
k2_binA = imcomplement(k2_binA); 
k2_binB = imcomplement(k2_binB);
% Create a Euclidean distance transform of the k2_bin
[distk2A,~] = bwdist(k2_binA,'euclidean');
[distk2B,~] = bwdist(k2_binB,'euclidean');
distk2A = double(distk2A)*0.6;
distk2B = double(distk2B)*0.6;
distk2 = (distk2A+distk2B)/2;
%% Rasterize images (compute average over every raster element
for i = 1:numel(imgName)
    [actinRaster,xActin,yActin] = rasterizeImage(Img(:,:,i),pixSize);
    meanIntensity = mean(actinRaster,[1 2]);
    Img_ras(:,:,i) = reshape(meanIntensity,yActin,xActin)'; % rasterized image, with "superpixels"
end
[gaussRaster,xGauss,yGauss] = rasterizeImage(gaussCurv,pixSize);
[meanRaster,xMean,yMean] = rasterizeImage(meanCurv,pixSize);
[k1Raster,xk1,yk1] = rasterizeImage(curvMap_k1,pixSize);
[k2Raster,xk2,yk2] = rasterizeImage(curvMap_k2,pixSize);
[distk2Raster,xDist,yDist] = rasterizeImage(distk2,pixSize);
meanGauss = mean(gaussRaster,[1 2]);
meanMean = mean(meanRaster,[1 2]);
meank1 = mean(k1Raster,[1 2]);
meank2 = mean(k2Raster,[1 2]);
meandistk2 = mean(distk2Raster,[1,2]);
gaussCurv_ras = reshape(meanGauss,yGauss,xGauss)';
gaussCurv_ras = repmat(gaussCurv_ras,1,1,length(imgName));
meanCurv_ras = reshape(meanMean,yMean,xMean)';
meanCurv_ras = repmat(meanCurv_ras,1,1,length(imgName));
k1_ras = reshape(meank1,yk1,xk1)';
k1_ras = repmat(k1_ras,1,1,length(imgName));
k2_ras = reshape(meank2,yk2,xk2)';
k2_ras = repmat(k2_ras,1,1,length(imgName));
distk2_ras = reshape(meandistk2,yDist,xDist)';
distk2_ras = repmat(distk2_ras,1,1,length(imgName));
%% Reshape data
Img_1D = reshape(Img_ras,[],1);
gaussCurv_1D = reshape(gaussCurv_ras,[],1);
meanCurv_1D = reshape(meanCurv_ras,[],1);
k1_1D = reshape(k1_ras,[],1);
k2_1D = reshape(k2_ras,[],1);
distk2_1D = reshape(distk2_ras,[],1);
%% Sign functions
sgnk1_1D = sign(k1_1D);
sgnk2_1D = sign(k2_1D);

k2zero_1D = Img_1D(sgnk2_1D==0);
k2neg_1D = Img_1D(sgnk2_1D<0);
k2pos_1D = Img_1D(sgnk2_1D>0);


k2zero_1D_2 = NaN*ones(size(Img_1D));
k2neg_1D_2 = NaN*ones(size(Img_1D));
k2pos_1D_2 = NaN*ones(size(Img_1D));

k2zero_1D_2(1:length(k2zero_1D)) = k2zero_1D;
k2neg_1D_2(1:length(k2neg_1D)) = k2neg_1D;
k2pos_1D_2(1:length(k2pos_1D)) = k2pos_1D;

%% Store all 1D arrays in a big array, specific to this caseExp
% k1, k2, K, H, distk2, k2neg, k2zero, k2pos, Img
Data_1D = [k1_1D,k2_1D,gaussCurv_1D,meanCurv_1D,distk2_1D,k2neg_1D_2,k2zero_1D_2,k2pos_1D_2,Img_1D];