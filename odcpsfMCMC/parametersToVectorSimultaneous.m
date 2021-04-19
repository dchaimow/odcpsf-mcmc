function v = parametersToVectorSimultaneous(...
    noise,rho,delta,epsilon,omega,theta,fwhm1,fwhm2)
v = [noise(:)' rho delta epsilon omega theta fwhm1 fwhm2];
end
