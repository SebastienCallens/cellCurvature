function ImgCropped = cropMasked(Img,mask)
% Remove black pixels around ROI that are the result from applying the mask
% to the original image
[r,c] = find(double(mask));
row1 = min(r);
row2 = max(r);
col1 = min(c);
col2 = max(c);
ImgCropped = Img(row1:row2, col1:col2);