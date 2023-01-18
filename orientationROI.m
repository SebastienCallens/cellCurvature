%%%% Detect orientation in ROI
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
substrate = 'Unduloid';
caseExp = 'D8ConvexDiff';

load('example_data/20191212_D8ConvexDiffBS4Unduloid_actin.mat');

% Normalize image 
Img = double(Img);
Img = Img/mean(Img,'all');

pixSize = 80; % Superpixel size for rasterization

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
% Now remove black pixels that result after applying the mask
ImgROI = cropMasked(Img,cropMask);
ImgROI_sharp = ImgROI;

%% Rasterize ROI and compute orientation
tic;
[actinRaster,xActin,yActin,xpos2,ypos2] = rasterizeImage(ImgROI_sharp,pixSize);
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
actinAngle = reshape(actinAngle,xActin,yActin)'; % Angle in degrees
actin_S = reshape(actin_S,xActin,yActin)';

%% Make plots
% Make two arrays with x and y positions of the quiver lines
xpos = 1:size(actinAngle,2);
xpos = repmat(xpos,size(actinAngle,1),1);
ypos = (size(actinAngle,1):-1:1)';
ypos = repmat(ypos,1,size(actinAngle,2));
% Make quiver plot with line length determined by orientation strength and
% actin intensity
figure(100)
subplot(1,2,1)
imagesc(ImgROI_sharp)
colormap(gray)
hold on
plot(reshape(xpos2,1,[]),reshape(ypos2,1,[]),'o','MarkerFaceColor',color3,'Color',color3,'MarkerSize',4)
axis equal
axis tight
subplot(1,2,2)
quiverAmplifier = 4;
quiver(xpos2,flipud(ypos2),quiverAmplifier*meanActin/max(max(meanActin)).*actin_S.*cosd(actinAngle),...
    quiverAmplifier*meanActin/max(max(meanActin)).*actin_S.*sind(actinAngle),'-','Color',color4,'LineWidth',1,'MaxHeadSize',0.3)
axis equal
axis tight
grid on

figure(200)
% Need to flip y pos to make arrows and figure match in the overlay plot
imagesc(ImgROI_sharp); colormap(gray);
hold on
quiver(xpos2,ypos2,quiverAmplifier*meanActin/max(max(meanActin)).*actin_S.*cosd(actinAngle),...
    -quiverAmplifier*meanActin/max(max(meanActin)).*actin_S.*sind(actinAngle),0.5,'-','Color',color3,'LineWidth',1,'MaxHeadSize',0.3)
for i = 1:length(reshape(xpos2,1,[]))
    xpatch = [xpos2(i)-pixSize/2,xpos2(i)+pixSize/2,xpos2(i)+pixSize/2,xpos2(i)-pixSize/2];
    ypatch = [ypos2(i)-pixSize/2,ypos2(i)-pixSize/2,ypos2(i)+pixSize/2,ypos2(i)+pixSize/2];
    patch('XData',xpatch,'YData',ypatch,'FaceColor','none','EdgeColor',color1,'LineWidth',1);
end
hold off
axis equal; set(gca,'visible','off')
