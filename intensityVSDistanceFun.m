function [] = intensityVSDistanceFun(substrate,caseExp,imgName)

color1 = [239,193,100]/255; % yellowish
color2 = [243,131,93]/255; % orangeish
color3 = [243,89,85]/255; % redish
color4 = [40,98,117]/255; % dark blueish
color5 = [0,67,76]/255; % darker blueish

% Load processed images
for i = 1:length(imgName)
    Imgtemp = load(strcat(imgName{i},'.mat'));
    Imgtemp = Imgtemp.Img;
    Imgtemp = double(Imgtemp);
    Imgtemp_norm = Imgtemp/mean(Imgtemp,'all');
    % Store non-normalized and normalized versions
    ImgIn(:,:,i) = Imgtemp;
    ImgIn_norm(:,:,i) = Imgtemp_norm;
end
Img = ImgIn_norm;
% Load curvature map
AA = load(strcat(substrate,'_k1Map.mat'));
BB = load(strcat(substrate,'_k2Map.mat'));
curvMap_k1 = AA.curvMap_k1;
curvMap_k2 = BB.curvMap_k2;
if strcmp(substrate,'Wavy')
    curvMap_k1 = flip(curvMap_k1,1);
    curvMap_k2 = flip(curvMap_k2,1);
end
%% Create Euclidean distance transform for the k2 curvature maps
k2_binA = imbinarize(curvMap_k2,-0.00000001);
k2_binB = k2_binA;
k2_binA(curvMap_k1==0) = 0; 
k2_binA = imcomplement(k2_binA);
k2_binB = imcomplement(k2_binB);
% Create a Euclidean distance transform of the k2_bin
[distk2A,~] = bwdist(k2_binA,'euclidean');
[distk2B,~] = bwdist(k2_binB,'euclidean');
distk2A = double(distk2A)*0.6; % convert to micron (projected)
distk2B = double(distk2B)*0.6;
distk2 = (distk2A+distk2B)/2;
distk2 = repmat(distk2,1,1,length(imgName));
%% Bin the image and the distance map
Img_1D = reshape(Img,[],1);
distk2_1D = reshape(distk2,[],1);
nbins = 30;
[distBinned,distEdges] = discretize(distk2_1D,nbins);
sampledData = [];
nSampling = 100;
for i=1:nbins
    meanIntDist(i) = median(Img_1D(distBinned==i));
    perc25Dist(i) = prctile(Img_1D(distBinned==i),25);
    perc75Dist(i) = prctile(Img_1D(distBinned==i),75);

    stdDist(i) = std(Img_1D(distBinned==i));
    nEl(i) = length(Img_1D(distBinned==i));
    tempSample = randsample(Img_1D(distBinned==i),nSampling);
    tempSampleDist = randsample(distk2_1D(distBinned==i),nSampling);
    sampledData = [sampledData;[tempSampleDist,tempSample]];

end

figure
locX = distEdges(1:end-1)+diff(distEdges(1:end))/2; % center points of bins
locX2 = distEdges(1:end-1); % left edges of bins
fill([locX2,fliplr(locX2)],[perc75Dist,fliplr(perc25Dist)],color5,'FaceAlpha',0.2,'edgecolor','none')
hold on
plot(locX2,meanIntDist,'-','Color',color5,'LineWidth',2)
xlabel('Projected distance from \kappa_2 \leq 0')
ylabel('Normalized intensity')

