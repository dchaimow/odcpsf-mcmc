function v = parametersToVector(...
    noise,rho,delta,epsilon,omega,theta,fwhm)
v = [noise(:)' rho delta epsilon omega theta fwhm];
end
