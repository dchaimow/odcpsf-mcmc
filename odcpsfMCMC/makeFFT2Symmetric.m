function As = makeFFT2Symmetric(A)
if isreal(A)
    As = (A + fft2DoubleFlip(A))/2;
else
    As = (A + fft2DoubleFlipConj(A))/2;
end
