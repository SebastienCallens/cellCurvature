%%%% Plot 2D heat map of intensity vs curvature
%%%% Sebastien Callens
clear; close all; clc
cmap = coolwarm(256);
cmap2 = fire(256);

caseExp1 = 'D8ConvexDiff';
caseExp2 = 'D8ConcaveDiff';
pixSize = 20;

saveFigQ=0; % Query whether to save figure or not

ConvexName = strcat('example_data/',caseExp1,'_all_',num2str(pixSize),'.xlsx');
ConcaveName = strcat('example_data/',caseExp2,'_all_',num2str(pixSize),'.xlsx');
ConvexData = xlsread(ConvexName);
ConcaveData = xlsread(ConcaveName);

k1Convex = ConvexData(:,1);
k2Convex = ConvexData(:,2);
KConvex = k1Convex.*k2Convex;
HConvex = 0.5*(k1Convex+k2Convex);
IConvex = ConvexData(:,9);
k1Concave = ConcaveData(:,1);
k2Concave = ConcaveData(:,2);
KConcave = k1Concave.*k2Concave;
HConcave = 0.5*(k1Concave+k2Concave);
IConcave = ConcaveData(:,9);

k1 = [k1Convex;k1Concave];
k2 = [k2Convex;k2Concave];
I = [IConvex;IConcave];
K = [KConvex;KConcave];
H = [HConvex;HConcave];

k1min = -5;
k1max = 20;
k2min = -k1max;
k2max = -k1min;

nBins = 40;
binEdgesk1 = linspace(k1min,k1max,nBins);
binEdgesk2 = linspace(k2min,k2max,nBins);
bink1 = discretize(k1,binEdgesk1);
bink2 = discretize(k2,binEdgesk2);

[X,Y] = meshgrid(binEdgesk1,binEdgesk2);

Idisc = NaN*ones(nBins);
for i = 1:nBins
    for j = 1:nBins
        logicalk1 = bink1==i;
        logicalk2 = bink2==j;
        logicalk1k2 = sum([logicalk1,logicalk2],2)==2;
        selectedI = I(logicalk1k2);
        Idisc(i,j) = median(selectedI);
    end
end

figure
colormap(cmap)
p = pcolor(X,Y,Idisc');
p.LineStyle = 'none';
hold on
plot([k1min,-k1min],[k1min,-k1min],'Color',[0.3,0.3,0.3],'LineStyle','--','LineWidth',1)
plot([0,k1max],[0,0],'Color',[0.3,0.3,0.3],'LineStyle','--','LineWidth',1)
plot([0,0],[k2min,0],'Color',[0.3,0.3,0.3],'LineStyle','--','LineWidth',1)
xlabel('$$\tilde{\kappa_1}$$ (-)','interpreter','LateX','FontSize',14);
ylabel('$$\tilde{\kappa_2}$$ (-)','interpreter','LateX','FontSize',14);
cc = colorbar;
caxis([0,2])
caxis([0,1.6]);
title(cc,'Median intensity')