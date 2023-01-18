function [rasterImg,xNum,yNum,xpos,ypos] = rasterizeImage(Img,pixSize)
%% Resize image depending on pixSize
xWidth = pixSize; % Number of pixels in x direction in one region
yWidth = pixSize; % Number of pixels in y direction in one region
xWidthTotal = size(Img,2);
yWidthTotal = size(Img,1);

% Total number of pixels might not be divisible by size of blocks. Select
% suitable starting points (and throw away some pixels at the edges)
xStart = ceil(mod(xWidthTotal,xWidth)/2);
if xStart == 0
    xStart = xStart+1;
end
yStart = ceil(mod(yWidthTotal,yWidth)/2);
if yStart == 0
    yStart = yStart+1;
end
xNum = floor(xWidthTotal/xWidth);
yNum = floor(yWidthTotal/yWidth);
% create resized image
Img_res = Img(yStart:yStart-1+yNum*yWidth,xStart:xStart-1+xNum*xWidth);

% reshape and permute
tempRaster = reshape(Img_res,yWidth,yNum,[]);
tempRaster_perm = permute(tempRaster,[1 3 2]);
rasterImg = reshape(tempRaster_perm,yWidth,xWidth,[]);

xpos = xStart+xWidth/2:xWidth:xStart+xNum*xWidth-xWidth/2;
xpos = repmat(xpos,yNum,1);

ypos = yStart+yWidth/2:yWidth:yStart+yNum*yWidth-yWidth/2;
ypos = repmat(ypos',1,xNum);
