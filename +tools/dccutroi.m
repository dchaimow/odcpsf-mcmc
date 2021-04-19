function [img,roi] = dccutroi(img,roi,d)
roi = roi>0.5;
if ~exist('d','var')
    d = 0;
end

[X,Y] = find(roi);
minX = max(1,min(X)-d);
minY = max(1,min(Y)-d);
maxX = min(size(img,1),max(X)+d);
maxY = min(size(img,2),max(Y)+d);

if ndims(img)>2
    sz = size(img);
    imgSemiFlat = reshape(img,[sz(1) sz(2) prod(sz(3:end))]);
    imgSemiFlat = imgSemiFlat(minX:maxX,minY:maxY,:);
    sz2 = size(imgSemiFlat);
    img = reshape(imgSemiFlat,[sz2(1) sz2(2) sz(3:end)]);
else
    img = img(minX:maxX,minY:maxY);
end
roi = roi(minX:maxX,minY:maxY);
end
