function B = mypadarray(A,padsize,direction)
if ndims(A)==2
   if exist('direction','var') && strcmp(direction,'post')
    origsize = size(A);
    newsize = origsize+padsize;
    B = zeros(newsize);
    B(1:origsize(1),1:origsize(2)) = A;
   else
    origsize = size(A);
    newsize = origsize+2*padsize;
    B = zeros(newsize);
    B(padsize(1)+(1:origsize(1)),padsize(2)+(1:origsize(2))) = A;      
   end
else
    error('direction needs to be post and ndims = 2');
end
