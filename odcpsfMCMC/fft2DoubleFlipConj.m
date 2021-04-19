function As = fft2DoubleFlipConj(A)
% generates the fft 'symmetric' version of a complex matrix A
% so that if ifft2(A) is real then As==A
As = conj(fft2DoubleFlip(A));