function Aup = dcupsample(A,factor)
fA = fft2(A);
fAp = ifftshift(mypadarray(fftshift(fA),(factor-1)*(size(A)/2)));
Aup = real(ifft2(fAp)*factor^2);
end
