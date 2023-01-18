%%%% Run many images for intensity vs distance
clear; close all; clc

substratelist = {'Unduloid'};
caseExp = 'D8ConvexDiff';

for  i = 1:length(substratelist)
    tic;
    substrate = substratelist{i};
    imgName = selectImgName(substrate,caseExp);
    intensityVSDistanceFun(substrate,caseExp,imgName);
    toc;
end