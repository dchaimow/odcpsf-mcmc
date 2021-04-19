function [noise,rho,delta,epsilon,omega,theta,fwhm] = ...
    vectorToParameters(v,sim,noMod)
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
fwhm    = v(:,nNoise+6);
end
