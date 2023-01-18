%%%% Intensity plot along the middle line of the rotationally symmetric
%%%% structures
%%%% Sebastien Callens
clear; close all; clc;
cmap = coolwarm(256);
cmap2 = fire(256);
color1 = [239,193,100]/255; % yellowish
color2 = [243,131,93]/255; % orangeish
color3 = [243,89,85]/255; % redish
color4 = [40,98,117]/255; % dark blueish
color5 = [0,67,76]/255; % darker blueish
%% Input

substrate = 'Unduloid';
caseImg = 'D8ConvexDiff';

switch substrate
    case 'Spheres'
        switch caseImg
            case 'D8ConvexDiff'
                ImgA = load('example_data/20191204_D8ConvexDiffBS2Spheres_actin.mat');
                ImgB = load('example_data/20191212_D8ConvexDiffS4BSpheres_actin.mat');
                ImgC = load('example_data/20200116_D8ConvexDiffS2Spheres_actin.mat');
                ImgD = load('example_data/20200121_D8ConvexDiffS4Spheres_actin.mat');
                ImgA = double(ImgA.Img);
                ImgB = double(ImgB.Img);
                ImgC = double(ImgC.Img);
                ImgD = double(ImgD.Img);
                Img = ImgA/mean(ImgA,'all')+ImgB/mean(ImgB,'all')+ImgC/mean(ImgC,'all')+ImgD/mean(ImgD,'all');
        end
    case 'Unduloid'
        switch caseImg
           
            case 'D8ConvexDiff'
                ImgA = load('example_data/20191204_D8ConvexDiffBS2Unduloid_actin.mat');
                ImgB = load('example_data/20191212_D8ConvexDiffBS4Unduloid_actin.mat');
                ImgC = load('example_data/20200116_D8ConvexDiffS2Unduloid_actin.mat');
                ImgD = load('example_data/20200121_D8ConvexDiffS4Unduloid_actin.mat');
                ImgA = double(ImgA.Img);
                ImgB = double(ImgB.Img);
                ImgC = double(ImgC.Img);
                ImgD = double(ImgD.Img);
                Img = ImgA/mean(ImgA,'all')+ImgB/mean(ImgB,'all')+ImgC/mean(ImgC,'all')+ImgD/mean(ImgD,'all');
            
        end
    case 'Catenoids'
        switch caseImg
            case 'D8ConvexDiff'
                ImgA = load('example_data/20191204_D8ConvexDiffBS2Catenoids_actin.mat');
                ImgB = load('example_data/20191212_D8ConvexDiffS4BCatenoids_actin.mat');
                ImgC = load('example_data/20200116_D8ConvexDiffS2Catenoid_actin.mat');
                ImgD = load('example_data/20200121_D8convexDiffS4Catenoid_actin.mat');
                ImgA = double(ImgA.Img);
                ImgB = double(ImgB.Img);
                ImgC = double(ImgC.Img);
                ImgD = double(ImgD.Img);
                Img = ImgA/mean(ImgA,'all')+ImgB/mean(ImgB,'all')+ImgC/mean(ImgC,'all')+ImgD/mean(ImgD,'all');
            
        end
    case 'Pseudospheres'
        switch caseImg
            case 'D8ConvexDiff'
                ImgA = load('example_data/20191204_D8ConvexDiffBS2Pseudospheres_actin.mat');
                ImgB = load('example_data/20191212_D8ConvexDiffS4BPseudospheres_actin.mat');
                ImgC = load('example_data/20200116_D8ConvexDiffS2Pseudospheres_actin.mat');
                ImgD = load('example_data/20200121_D8ConvexDiffS4Pseudospheres_actin.mat');
                ImgA = double(ImgA.Img);
                ImgB = double(ImgB.Img);
                ImgC = double(ImgC.Img);
                ImgD = double(ImgD.Img);
                Img = ImgA/mean(ImgA,'all')+ImgB/mean(ImgB,'all')+ImgC/mean(ImgC,'all')+ImgD/mean(ImgD,'all');
           
        end
    case 'Cylinder'
        switch caseImg
            case 'D8ConvexDiff'
                ImgA = load('example_data/20191204_D8ConvexDiffBS2Cylinder_actin.mat');
                ImgB = load('example_data/20191212_D8ConvexDiffS4BCylinder_actin.mat');
                ImgC = load('example_data/20200121_D8ConvexDiffS4Cylinder_actin.mat');
                ImgA = double(ImgA.Img);
                ImgB = double(ImgB.Img);
                ImgC = double(ImgC.Img);
                Img = ImgA/mean(ImgA,'all')+ImgB/mean(ImgB,'all')+ImgC/mean(ImgC,'all');
            
        end
end
%% Load curvature map and color image
load(strcat(substrate,'_k1Map.mat'));
load(strcat(substrate,'_k2Map.mat'));
if strcmp(substrate,'Wavy')
    curvMap_k1 = flip(curvMap_k1,1);
    curvMap_k2 = flip(curvMap_k2,1);
end
% Normalize curvature by radius of sphere 
curvMap_k1 = curvMap_k1*180;
curvMap_k2 = curvMap_k2*180;

%% Extract intensity along middle line
% Determine top and bottom point for every structure
top = 1;
bottom = size(Img,1);
% Get midline (column number) of the structure
midLine = round(size(Img,2)/2);
% Get intensity vector along midline, from top to bottom
midInt = Img(top:bottom,midLine);
midk1 = curvMap_k1(top:bottom,midLine);
midk2 = curvMap_k2(top:bottom,midLine);
%% Plot
figure
subplot(2,1,1)
plot(0.6*(1:1:length(midInt)),midInt,'Color',color1)
hold on
plot(0.6*(1:1:length(midInt)),smooth(midInt,80),'Color',color3,'LineWidth',2)
ylabel('Normalized intensity (a.u.)')
set(gca,'xtick',[])
axis tight
subplot(2,1,2)
plot(0.6*(1:1:length(midInt)),midk1,'Color',color4,'LineWidth',2)
hold on
plot(0.6*(1:1:length(midInt)),midk2,'Color',color2,'LineWidth',2)
xlabel('Position (\mum)')
ylabel('Principal curvature (-)')
legend('\kappa_1','\kappa_2')
legend('boxoff')
axis tight


