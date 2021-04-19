function As = fft2DoubleFlip(A)
% generates the fft 'symmetric' version of a real matrix A
% so that if ifft2(A) is real then As==A
M = size(A,1);
N = size(A,2);
As = A(mod(M:-1:1, M) + 1, mod(N:-1:1, N) + 1);