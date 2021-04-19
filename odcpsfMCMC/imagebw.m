function imagebw(img,roi,range)
if ~exist('range','var')
    range = [min(img(:)) max(img(:))];
end
d = range(2)-range(1);
img = (img-range(1))/d;
img(roi==false) = 1;
img(img<0) = 0;
img(img>1) = 1;
image(repmat(img,[1 1 3]));
end
