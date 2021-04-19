function data = dcunflat(flatdata,roi)
roi = roi>0.5;
nT = size(flatdata,1);
data = zeros([numel(roi) nT]);
data(roi==1,:) = flatdata';
data = reshape(data,[size(roi) nT]);
end
