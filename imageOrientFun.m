function [maxAngle,orientStrength,P] = imageOrientFun(Img)
%% Perform windowing of the image to avoid artefacts from assumed periodicity
window1 = cos(linspace(-pi/2,pi/2,size(Img,1)));
window2 = cos(linspace(-pi/2,pi/2,size(Img,2)));
window = window1'*window2;

Img = im2single(Img).*window;
%% Calculate Fast Fourier Transform (FFT) Based on: Fourier analysis and automated measurement of cell and fiber angular orientation distributions
Y = fft2(Img);
P = log(abs(fftshift(Y)) + 1); % The term +1 in the logarithm is to account for issues with 0 (logarithm is just to make the graphs more "readable")
%% Calculate orientation with P
% We calculate the principal orientation by rotating the image through all
% angles and summing the power values of a central band (a range of
% columns). The angle where this becomes maximal is the orientation angle. 
theta = linspace(0,360,361); % We do the orientation over 360 degrees. This is to make sure we can also detect potential peaks at the endpoints
Psum = zeros(length(theta),1);
for i = 1:length(theta)
    Imgrot = imrotate(P,theta(i));
    Imgrot(Imgrot==0) = NaN; % Set the zeros that result after rotating the image to NaN (to exclude those pixels)
    
    midCol = round(size(Imgrot,2)/2); % Central column of the data
    colRange = 1; % Range of columns around central column that we will consider
    Imgrot_slice = Imgrot(:,midCol-colRange:midCol+colRange);
    Psum(i) = nansum(nansum(Imgrot_slice));
    numP = sum(sum(~isnan(Imgrot_slice))); % Number of elements in the slice that are not NaNs
    Psum(i) = Psum(i)/numP;
    
end
Psum = Psum/sum(Psum);
% Use peak analysis to find principal orientation and strength of
% orientation
[P_peaks,P_locs,P_width,P_prominence] = findpeaks(Psum,theta,'MinPeakDistance',5,'Annotate','extents'); % We ignore the first five datapoints because we extended the data to 185 degrees to find the "end" peak
firstHalf = P_prominence;
diffWithMax = firstHalf-max(firstHalf);
[~,minIdx] = min(abs(diffWithMax));
maxAngle = 180-P_locs(firstHalf==max(firstHalf));

if maxAngle<0
    maxAngle = maxAngle+180;
elseif maxAngle>180
    maxAngle = maxAngle-180;
end

if length(maxAngle)~=1
    maxAngle = max(abs(maxAngle));
end
if isempty(maxAngle)
    maxAngle = NaN;
end
orientStrength = max(P_prominence)/(sum(P_prominence)/2); 
if isempty(orientStrength)
   orientStrength = 0;
end
