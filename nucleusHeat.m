%%%% Generate heat maps of nuclei
%%%% Sebastien Callens
clear; close all; clc

clear; close all; clc;
cmap = coolwarm(256);
cmap2 = fire(256);

color1 = [239,193,100]/255; % yellowish
color2 = [243,131,93]/255; % orangeish
color3 = [243,89,85]/255; % redish
color4 = [40,98,117]/255; % dark blueish
color5 = [0,67,76]/255; % darker blueish
%% Input
%1. Specify type of substrate
%2. Specify case (e.g. D8Diff)
%3. Specify which images (e.g. ImgA-ImgD)
%4. Run script

substrate = 'Unduloid';
caseImg = 'D5ConvexDiff';
saveFigQ = 0; % query to save figure or not

nucleiFolder = 'example_data\';

fileList = {'20191209_D5ConvexDiffS2Unduloid_nucleus_data',...
    '20200203_D5ConvexDiffS2Unduloid_nucleus_data',...
    '20200204_D5ConvexDiffS4Unduloid_nucleus_data'};

nuclei = [];
% Combine nuclei centroids
for i = 1:length(fileList)
    data = xlsread(strcat(nucleiFolder,fileList{i},'.xlsx'));
    nucleiAll{i} = data(:,3:4);
    nuclei = [nuclei;data(:,3:4)];
end
% Plot nuclei centroids
figure;
plot(nuclei(:,1),nuclei(:,2),'r.')
axis equal; axis tight;


% Round nuclei centroid positions for matching with curvature and distance
% maps
nucleiRound = round(nuclei);

%% Distance map
load(strcat(substrate,'_k1Map.mat'));
load(strcat(substrate,'_k2Map.mat'));
if strcmp(substrate,'Wavy')
    curvMap_k1 = flip(curvMap_k1,1);
    curvMap_k2 = flip(curvMap_k2,1);
end
% Normalize curvature by radius of sphere (we have to pick some length
% measure to normalize, so sphere radius is picked as reference)
curvMap_k1 = curvMap_k1/180;
curvMap_k2 = curvMap_k2/180;

% Binarize the k2 curvature map with a threshold of slightly below 0 (all values larger
% than 0 will be set equal to 1).
% We create two maps: one for the distance from k2 <= 0 and one for the
% distance from k2 <0. We then take the average of those two maps. The
% reason is to include some "attraction effect" from the negative k2
% regions, while still making sure that k2=0 is not as bad as k2>0
k2_binA = imbinarize(curvMap_k2,-0.00000001);
k2_binB = k2_binA;
k2_binA(curvMap_k1==0) = 0; % This is to also include the flat points. I just want to make sure that the points where k1>0 and k2=0 are not included
k2_binA = imcomplement(k2_binA); % I want the values of 0 or smaller to be 1, and the other ones 0
k2_binB = imcomplement(k2_binB);
% Create a Euclidean distance transform of the k2_bin
[distk2A,~] = bwdist(k2_binA,'euclidean');
[distk2B,~] = bwdist(k2_binB,'euclidean');
%distk2A = distk2A/max(max(distk2A)); % normalize w.r.t. max value
%distk2B = distk2B/max(max(distk2B));
% convert to double (just to make sure it is the same format as Img)
distk2A = double(distk2A)*0.6;
distk2B = double(distk2B)*0.6;
distk2 = (distk2A+distk2B)/2;
% distk2 = distk2A;

%% Map nuclei to curvature value and distk2 value
for i = 1:size(nucleiRound,1)
    nucleiData(i,:) = [nucleiRound(i,1),nucleiRound(i,2),curvMap_k1(nucleiRound(i,2),nucleiRound(i,1)),curvMap_k2(nucleiRound(i,2),nucleiRound(i,1)),distk2(nucleiRound(i,2),nucleiRound(i,1))];
end
nBins = 20;
[distk2_counts,distk2_edges] = histcounts(distk2,nBins);
nuclei_counts = histcounts(nucleiData(:,5),nBins);
nuclei_density = nuclei_counts./distk2_counts;
distk2_mids = 0.5*(distk2_edges(1:end-1)+distk2_edges(2:end));

figure;
% plot(distk2_mids,nuclei_density,'o','Color',color5,'MarkerFaceColor',color5,'LineWidth',2)
f = fit(distk2_mids',nuclei_density','smoothingspline','Exclude',numel(nuclei_density));
fig_f = plot(f,distk2_mids',nuclei_density');
set(fig_f,'Color',color5,'LineWidth',2,'MarkerSize',30)
set(gca,'Fontsize',16)
legend('location','SouthWest')
xlabel('Projected distance from \kappa_2 \leq 0')
ylabel('Nuclei density')


%% Split and collapse nuclei data for repetitive units
switch substrate
    case 'Unduloid'
        top = 86; %Row number of top cut-line
        bottom = 4309; %Row number of bottom cut-line
        numRep = 5; %Number of repetitions
        lengthRep = round((bottom-top)/numRep);
        plot_row = 2;
        plot_col = 3;
    case 'Spheres'
%         top = 100;
%         bottom = 4122;
        top = 65;
        bottom = 4115;
        numRep = 7;
        lengthRep = round((bottom-top)/numRep);
        plot_row = 2;
        plot_col = 4;
    case 'Catenoids'
%         top = 74;
%         bottom = 3878;
        top = 90;
        bottom = 3870;
        numRep = 15/2;
        lengthRep = round((bottom-top)/numRep);
        plot_row = 4;
        plot_col = 4;
    case 'Pseudospheres'
        top = 86;
        bottom = 3961;
        numRep = 7;
        lengthRep = round((bottom-top)/numRep);
        plot_row = 2;
        plot_col = 4;
    case 'Cylinder'
        top = 78;
        bottom = 4302;
        numRep = 5;
        lengthRep = round((bottom-top)/numRep);
        plot_row = 2;
        plot_col = 3;
    case 'Wavy'
        top = 255;
        bottom = 4214;
        numRep = 2;
        lengthRep = round((bottom-top)/numRep);
        plot_row = 1;
        plot_col = 3;
        
end

% Remove top and bottom
nucleiCrop = nuclei(nuclei(:,2)>top&nuclei(:,2)<bottom,:);
% Collapse units, using modulo. Divisor is distance between max and min,
% divided by numRep
if strcmp(substrate,'Wavy')
    nucleiRep = [nucleiCrop(:,1),mod(nucleiCrop(:,2),round((max(nucleiCrop(:,2))-min(nucleiCrop(:,2)))/numRep))];
else
    nucleiRep = [nucleiCrop(:,1),mod(nucleiCrop(:,2),round((max(nuclei(:,2))-min(nuclei(:,2)))/numRep))];
end

% Determine xWidth and yWidth of image
xWidth = max(nucleiRep(:,1))-min(nucleiRep(:,1));
yWidth = max(nucleiRep(:,2))-min(nucleiRep(:,2));
AR = yWidth/xWidth; % aspect ratio, necessary for ybin vs xbin (to get square bins)
xBins = 100;
yBins = round(AR*xBins);
% Create 2D histogram count
N = histcounts2(nucleiRep(:,1),nucleiRep(:,2),[xBins,yBins]);
N = N'; % for plotting

% Create Gaussian filter matrix
[xG, yG] = meshgrid(-3:3);
sigma = 1.5;
g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
g = g./sum(g(:));

% Convolve Gaussian filter with the N array and plot
figure
subplot(1,2,1)
imagesc(conv2(N,g,'same')); colormap(cmap2); axis equal; axis tight; colorbar
% Create version of nuclei centroid coordinates
scaleFactor = (max(nucleiRep(:,1))-min(nucleiRep(:,1)))/xBins;
nucleiScaled = nucleiRep/scaleFactor;
subplot(1,2,2)
% imagesc(conv2(N,g,'same')); colormap(cmap2); axis equal; axis tight; colorbar
% hold on
% plot(nucleiScaled(:,1),nucleiScaled(:,2),'w.'); 
plot(nucleiRep(:,1),nucleiRep(:,2),'k.'); axis equal; axis tight;

figure
imagesc(conv2(N,g,'same')); colormap(cmap2); axis equal; axis tight; colorbar
if strcmp(substrate,'Spheres') && strcmp(caseImg,'D8ConvexDiff')
    caxis([0,3])
end
axis off
saveStr = strcat('Figures/nucleusHeat_',substrate,'_',caseImg);
if saveFigQ==1
    print(saveStr,'-dsvg');
end