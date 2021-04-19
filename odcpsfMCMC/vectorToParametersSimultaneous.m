function [noise,rho,delta,epsilon,omega,theta,fwhm1,fwhm2] = ...
    vectorToParametersSimultaneous(v,sim,noMod)
nNoise  = sim.N1 * sim.N2;
noise   = squeeze(reshape(v(:,1:nNoise),...
    [size(v,1) sim.N1 sim.N2]));
rho     = v(:,nNoise+1);
delta   = v(:,nNoise+2);
epsilon = v(:,nNoise+3);
omega   = v(:,nNoise+4);
if exist('noMod','var') && noMod
    theta   = v(:,nNoise+5);
else
    theta   = mod(v(:,nNoise+5),pi);
end
fwhm1   = v(:,nNoise+6);
fwhm2   = v(:,nNoise+7);
end
