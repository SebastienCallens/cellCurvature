%%%% Frequency map generation
%%%% Sebastien Callens
clear; close all; clc;
cmap2 = fire(256);
%% Input
%1. Specify type of substrate
%2. Specify case (e.g. D8Diff)
%3. Specify images
%4. Run script

substrate = 'Unduloid';
caseImg = 'D8ConvexDiff';
saveFigQ = 0; % query to save figure or not
imgName = {'example_data/20191204_D8ConvexDiffBS2Unduloid_actin',...
           'example_data/20191212_D8ConvexDiffBS4Unduloid_actin',...
           'example_data/20200116_D8ConvexDiffS2Unduloid_actin',...
           'example_data/20200121_D8ConvexDiffS4Unduloid_actin'};

for i = 1:length(imgName)
    Imgtemp = load(strcat(imgName{i},'.mat'));
    Imgtemp = Imgtemp.Img;
    Imgtemp = double(Imgtemp);
    Imgtemp_norm = Imgtemp/mean(Imgtemp,'all');
    ImgAll(:,:,i) = Imgtemp_norm;
end
Img = sum(ImgAll,3);
Img = Img/mean(Img,'all');

%% Split images
switch substrate
    case 'Unduloid'
        top = 86; 
        bottom = 4309; 
        numRep = 5; 
        lengthRep = round((bottom-top)/numRep);
    case 'Spheres'
        top = 65;
        bottom = 4115;
        numRep = 7;
        lengthRep = round((bottom-top)/numRep);
    case 'Catenoids'
        top = 90;
        bottom = 3870;
        numRep = 15/2;
        lengthRep = round((bottom-top)/numRep);
    case 'Pseudospheres'
        top = 86;
        bottom = 3961;
        numRep = 7;
        lengthRep = round((bottom-top)/numRep);
    case 'Cylinder'
        top = 78;
        bottom = 4302;
        numRep = 5;
        lengthRep = round((bottom-top)/numRep);
    case 'Wavy'
        top = 255;
        bottom = 4214;
        numRep = 2;
        lengthRep = round((bottom-top)/numRep);
end

for i = 1:numRep
    splitImg(:,:,i) = Img(top+(i-1)*lengthRep:top+i*lengthRep,:);
end

freqMap = sum(splitImg,3);

%% Plotting
freqFig = figure;
colormap(cmap2)
imagesc(freqMap)
axis equal
axis tight
set(gca,'visible','off')



