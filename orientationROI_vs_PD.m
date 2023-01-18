%%%% Detect orientation in ROI and compare to pd
%%%% Sebastien Callens
clear; close all; clc
cmap = coolwarm(256);
color1 = [239,193,100]/255; % yellowish
color2 = [243,131,93]/255; % orangeish
color3 = [243,89,85]/255; % redish
color4 = [40,98,117]/255; % dark blueish
color5 = [0,67,76]/255; % darker blueish
tic;
%% Input
% Load image
substrate = 'Unduloid';
load('example_data/20191212_D8ConvexDiffBS4Unduloid_actin.mat');
load(strcat(substrate,'_p1x.mat'));
load(strcat(substrate,'_p1y.mat'));
load(strcat(substrate,'_p2x.mat'));
load(strcat(substrate,'_p2y.mat'));
% Normalize image
Img = double(Img);
Img = Img/mean(Img,'all');

switch substrate
    case 'Unduloid'
        p1x = curvMap_p1y;
        p1y = -curvMap_p1x;
        p2x = curvMap_p2y;
        p2y = curvMap_p2x;
    case 'Pseudospheres'
        p1x = curvMap_p1x;
        p1y = -curvMap_p1y;
        p2x = curvMap_p2x;
        p2y = curvMap_p2y;
    case 'Cylinder'
        p1x = curvMap_p1y;
        p1y = curvMap_p1x;
        p2x = curvMap_p2y;
        p2y = curvMap_p2x;
    case 'Catenoids'
        p1x = -curvMap_p2x;
        p1y = curvMap_p2y;
        p2x = curvMap_p1x;
        p2y = curvMap_p2y;
    case 'Wavy'
        p1x = curvMap_p1x;
        p1y = curvMap_p1y;
        p2x = curvMap_p2x;
        p2y = curvMap_p2y;
end


pixSize = 80; % Superpixel size for rasterization
actinTol = 0.3; % tolerance (0-1) for actin intensity to be included in comparison to pd

toc;
%% Create ROI and masked image
figure
imagesc(Img); axis equal; axis tight;
rectROI = drawrectangle('Rotatable',true,'LineWidth',1);
queryRect = input('Press enter to continue...'); % allows us to adjust the rectangle before continuing

xLT = rectROI.Position(1); % x-coordinate of top left point of the non-rotated rectangle
yLT = rectROI.Position(2);
xRT = xLT+rectROI.Position(3);
yRT = yLT;
xRB = xRT;
yRB = yLT+rectROI.Position(4);
xLB = xLT;
yLB = yRB;
xCoord = [xLT;xRT;xRB;xLB];
yCoord = [yLT;yRT;yRB;yLB];

% Mask images
cropMask = roipoly(Img,xCoord',yCoord');
ImgROI = cropMasked(Img,cropMask);

p1xROI = cropMasked(p1x,cropMask);
p1yROI = cropMasked(p1y,cropMask);
p2xROI = cropMasked(p2x,cropMask);
p2yROI = cropMasked(p2y,cropMask);

% Sharpen image for better detection of stress fibers
ImgROI_sharp = imsharpen(ImgROI);

%% Rasterize ROI and compute orientation
tic;
[actinRaster,xActin,yActin,xpos2,ypos2] = rasterizeImage(ImgROI_sharp,pixSize);
[p1xRaster,x_p1x,y_p1x] = rasterizeImage(p1xROI,pixSize);
[p1yRaster,x_p1y,y_p1y] = rasterizeImage(p1yROI,pixSize);
[p2xRaster,x_p2x,y_p2x] = rasterizeImage(p2xROI,pixSize);
[p2yRaster,x_p2y,y_p2y] = rasterizeImage(p2yROI,pixSize);
toc;
actinAngle = zeros(1,size(actinRaster,3));
actin_S = zeros(size(actinAngle));
for i = 1:size(actinRaster,3)
    [actinAngle(i),actin_S(i),Pplot(:,:,i)] = imageOrientFun(actinRaster(:,:,i));
end
toc;

%% Compute means per raster element
meanActin = mean(actinRaster,[1 2]);
meanActin = reshape(meanActin,xActin,yActin)';
meanp1x = mean(p1xRaster,[1,2]);
meanp1x = reshape(meanp1x,x_p1x,y_p1x)';
meanp1y = mean(p1yRaster,[1,2]);
meanp1y = reshape(meanp1y,x_p1y,y_p1y)';
meanp2x = mean(p2xRaster,[1,2]);
meanp2x = reshape(meanp2x,x_p2x,y_p2x)';
meanp2y = mean(p2yRaster,[1,2]);
meanp2y = reshape(meanp2y,x_p2y,y_p2y)';
actinAngle = reshape(actinAngle,xActin,yActin)'; % Angle in degrees
actin_S = reshape(actin_S,xActin,yActin)';
% Set strength to zero if value is below a threshold
actin_S(actin_S<0.3)=0;

%% Compare SF angle with PD
p1Angle = rad2deg(angle(meanp1x+1i*meanp1y));
p2Angle = rad2deg(angle(meanp2x+1i*meanp2y));
% make sure all angles are in the positive range: 0-180
p1Angle(p1Angle<0)=180+p1Angle(p1Angle<0);
p2Angle(p2Angle<0)=180+p2Angle(p2Angle<0);
actinAngle(actinAngle<0)=180+actinAngle(actinAngle<0);
% Calculate differences between SF and PD
diff_p1Actin = abs(p1Angle-actinAngle);
diff_p2Actin = abs(p2Angle-actinAngle);
% % If difference is larger than 90, take complement
diff_p1Actin(diff_p1Actin>90) = 180-diff_p1Actin(diff_p1Actin>90);
diff_p2Actin(diff_p2Actin>90) = 180-diff_p2Actin(diff_p2Actin>90);
% Exclude the points outside the substrate (points where both meanp1x and
% meanp1y are zero (or where meanp2x and meanp2y are zero))
diff_p1Actin(abs(meanp1x)+abs(meanp1y)<0.0001)=NaN;
diff_p2Actin(abs(meanp1x)+abs(meanp1y)<0.0001)=NaN;
% Exclude the points where the strength or intensity is too low
diff_p1Actin(actin_S==0)=NaN;
diff_p2Actin(actin_S==0)=NaN;
diff_p1Actin(meanActin<actinTol*mean(meanActin,'all'))=NaN;
diff_p2Actin(meanActin<actinTol*mean(meanActin,'all'))=NaN;

% Calculate degree of alignment: DA = 1-diff/90 (0 if perpendicular, 1 if
% aligned)
DA_p1 = 1-diff_p1Actin/90;
DA_p2 = 1-diff_p2Actin/90;

% Create vectors with p1Angle, p2Angle and actinAngle of the "valid" points
p1Anglevec = p1Angle(~isnan(diff_p1Actin));
p1Anglevec = p1Anglevec(:);
p2Anglevec = p2Angle(~isnan(diff_p1Actin));
p2Anglevec = p2Anglevec(:);
actinAnglevec = actinAngle(~isnan(diff_p1Actin));
actinAnglevec = actinAnglevec(:);

%% Make plots
% Make two arrays with x and y positions of the quiver lines
xpos = 1:size(actinAngle,2);
xpos = repmat(xpos,size(actinAngle,1),1);
ypos = (size(actinAngle,1):-1:1)';
ypos = repmat(ypos,1,size(actinAngle,2));
% Make quiver plot with line length determined by orientation strength and
% actin intensity
figure
subplot(2,2,1)
imagesc(flipud(ImgROI_sharp))
axis equal
axis tight
subplot(2,2,2)
quiverAmplifier = 4;
quiver(xpos2,ypos2,-quiverAmplifier*meanActin/max(max(meanActin)).*actin_S.*cosd(actinAngle),...
    quiverAmplifier*meanActin/max(max(meanActin)).*actin_S.*sind(actinAngle),'-','Color',color4,'LineWidth',1,'MaxHeadSize',0.3)
axis equal
axis tight
grid on
subplot(2,2,3)
quiver(xpos2,ypos2,-meanp1x,meanp1y,'-','Color',color3,'LineWidth',1,'MaxHeadSize',0.3)
axis equal
axis tight
grid on
subplot(2,2,4)
quiver(xpos2,ypos2,-meanp2x,meanp2y,'-','Color',color3,'LineWidth',1,'MaxHeadSize',0.3)
axis equal
axis tight
grid on

% plot angle difference histogram
figure
binEdges = linspace(0,1,20);
histogram(DA_p1,binEdges,'Normalization','pdf','EdgeColor','k','FaceColor',color5)
xlabel('DA')
ylabel('Probability density')
xlim([0,1])