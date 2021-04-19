function [odc,mri,E,...
    dE_dNoiseij,dE_dRho,dE_dDelta,dE_dEpsilon,...
    dE_dTheta,dE_dOmega,dE_dFwhm] = computeODCModel(sim,...
    noise,rho,delta,epsilon,theta,omega,beta,fwhm,...
    data,sigmaData,roi,filterCutoffs,limitVector)
% odc = computeODCModel(sim,...
%    noise,rho,delta,epsilon,theta,omega):
% simulates an pattern
%
% [odc,mri] = computeODCModel(sim,...
%    noise,rho,delta,epsilon,theta,omega,beta,fwhm):
% simulates an odc and an mri pattern
%
% [odc,mri,E] = computeODCModel(sim,...
%    noise,rho,delta,epsilon,theta,omega,beta,fwhm,...
%    data,sigmaData,roi):
% returns in addition the energy of the model given data
%
% [odc,mri,E,...
%    dE_dNoiseij,dE_dRho,dE_dDelta,dE_dEpsilon,...
%    dE_dTheta,dE_dOmega,dE_dFwhm] = computeODCModel(sim,...
%    noise,rho,delta,epsilon,theta,omega,beta,fwhm,...
%    data,sigmaData,roi):
% returns in addition the gradient of the energy

% limitVector: limits parameters (for now only omega and rho)
%{rhoLimit, omegaLimit}
% parameterLimit = [lowerLimit upperLimit]
% no lower limit: lowerLimit = -Inf
% no upper limit: upperLimit = Inf
if exist('limitVector','var')
    rhoLimit   = limitVector{1};
    omegaLimit = limitVector{2};
else
    rhoLimit   = [-Inf Inf];
    omegaLimit = [-Inf Inf];
end


% set all prior parameters
rho_mean        = 0.57;
rho_sigma       = 0.10;
delta_mean      = 0.24;
delta_sigma     = 0.02;
epsilon_mean    = 0.67;
epsilon_sigma   = 0.19;
noise_sigma     = 0.1;  % this is also being used for filter normalization

%% calculate odc filter shape
rMinusRho = sim.kr-rho;
rPlusRho  = sim.kr+rho;

phiMinusTheta = sim.kphi-theta;
cosPhiMinusTheta = cos(phiMinusTheta);
sinPhiMinusTheta = sin(phiMinusTheta);

FODCNotNormalizedRPlus  = exp(-rMinusRho.^2/(2*delta^2));
FODCNotNormalizedRMinus = exp(-rPlusRho.^2/(2*delta^2));

FODCNotNormalizedAng1   = exp( cosPhiMinusTheta/epsilon^2);
FODCNotNormalizedAng2   = exp(-cosPhiMinusTheta/epsilon^2);

FODCNotNormalizedR   = FODCNotNormalizedRPlus + FODCNotNormalizedRMinus;
FODCNotNormalizedAng = FODCNotNormalizedAng1  + FODCNotNormalizedAng2;

FODCNotNormalized = FODCNotNormalizedR.*FODCNotNormalizedAng;

% make sure filter is symmetric:
% (because of numerical inaccuracies in sim.kphi)
FODCNotNormalized = makeFFT2Symmetric(FODCNotNormalized);

C = sqrt(mean(FODCNotNormalized(:).^2)) * noise_sigma;
FODC = FODCNotNormalized/C;

%% compute odc part
noiseF = fft2(noise);
odcLin = ifft2(FODC.*noiseF);

if omega==0
    odc = sign(odcLin);
elseif omega==Inf
    odc = odcLin;
else
    s_hat  = 1./(1+exp(-odcLin/omega)); % sigmoid(odcLin/omega)
    odc    = 2 * (s_hat - 0.5);
end

if nargout == 1
    return;
end

%% calculate BOLD and MRI filter shape
sigmaPSF = fwhm/(2*sqrt(2*log(2)));
FMTF = beta * exp(-2*pi^2*sigmaPSF^2*sim.mrikrsq)/sim.upSampleFactor^2;

if exist('filterCutoffs','var') && ~isempty(filterCutoffs)
    FMTF((sim.mrikrsq<max(filterCutoffs(1),0)^2)|...
        (sim.mrikrsq>filterCutoffs(2)^2)) = 0;
end

%% compute BOLD and MRI part
odcF   = fft2(odc);
mri    = ifft2(FMTF.*...
    makeFFT2Symmetric(odcF(sim.MRIIndices1,sim.MRIIndices2)));

if nargout == 2
    return;
end

%% compute E
inROI = (roi(:) == 1);
dRho    = rho - rho_mean;
dDelta  = delta - delta_mean;
dEpsilon= epsilon - epsilon_mean;
d       = (mri - data).*roi;
if fwhm<0 || omega<0 || delta<0 || epsilon<0 || rho<0 || ...
        rho<rhoLimit(1) || rho>rhoLimit(2) || ...
        omega<omegaLimit(1) || omega>omegaLimit(2);
    E = Inf;
else
    rhoE    = dRho^2/(2*rho_sigma^2);
    deltaE  = dDelta^2/(2*delta_sigma^2);
    epsilonE= dEpsilon^2/(2*epsilon_sigma^2);
    noiseE  = sum(noise(:).^2)./(2*noise_sigma^2);
    dataE   = sum(d(inROI).^2)./(2*sigmaData^2);
    E = dataE + noiseE + rhoE + deltaE + epsilonE;
end

if nargout == 3
    return;
end

%% odc filter derivatives (to be used for energy gradients)
dFODCNotNormalizedRPlus_dRho    = ...
    FODCNotNormalizedRPlus.*rMinusRho/delta^2;
dFODCNotNormalizedRMinus_dRho   = ...
    -FODCNotNormalizedRMinus.*rPlusRho/delta^2;
dFODCNotNormalizedRPlus_dDelta  = ...
    dFODCNotNormalizedRPlus_dRho.*rMinusRho/delta;
dFODCNotNormalizedRMinus_dDelta = ...
    -dFODCNotNormalizedRMinus_dRho.*rPlusRho/delta;

dFODCNotNormalizedAng1_dTheta = ...
    FODCNotNormalizedAng1.*sinPhiMinusTheta/epsilon^2;
dFODCNotNormalizedAng2_dTheta = ...
    FODCNotNormalizedAng2.*-sinPhiMinusTheta/epsilon^2;

dFODCNotNormalizedAng1_dEpsilon = ...
    FODCNotNormalizedAng1.*-2.*cosPhiMinusTheta/epsilon^3;
dFODCNotNormalizedAng2_dEpsilon = ...
    FODCNotNormalizedAng2.*2.*cosPhiMinusTheta/epsilon^3;

dFODCNotNormalizedR_dRho        = ...
    dFODCNotNormalizedRPlus_dRho+dFODCNotNormalizedRMinus_dRho;
dFODCNotNormalizedR_dDelta      = ...
    dFODCNotNormalizedRPlus_dDelta+dFODCNotNormalizedRMinus_dDelta;
dFODCNotNormalizedAng_dTheta    = ...
    dFODCNotNormalizedAng1_dTheta+dFODCNotNormalizedAng2_dTheta;
dFODCNotNormalizedAng_dEpsilon  = ...
    dFODCNotNormalizedAng1_dEpsilon+dFODCNotNormalizedAng2_dEpsilon;

dFODCNotNormalized_dRho     = ...
    dFODCNotNormalizedR_dRho.*FODCNotNormalizedAng;
dFODCNotNormalized_dDelta   = ...
    dFODCNotNormalizedR_dDelta.*FODCNotNormalizedAng;
dFODCNotNormalized_dTheta   = ...
    FODCNotNormalizedR.*dFODCNotNormalizedAng_dTheta;
dFODCNotNormalized_dEpsilon = ...
    FODCNotNormalizedR.*dFODCNotNormalizedAng_dEpsilon;

dC_dRho = ...
    noise_sigma^2*sum(FODCNotNormalized(:).*dFODCNotNormalized_dRho(:))...
    /(C*sim.N1*sim.N2);
dC_dDelta = ...
    noise_sigma^2*sum(FODCNotNormalized(:).*dFODCNotNormalized_dDelta(:))...
    /(C*sim.N1*sim.N2);
dC_dTheta = ...
    noise_sigma^2*sum(FODCNotNormalized(:).*dFODCNotNormalized_dTheta(:))...
    /(C*sim.N1*sim.N2);
dC_dEpsilon = ...
    noise_sigma^2*sum(FODCNotNormalized(:).*dFODCNotNormalized_dEpsilon(:))...
    /(C*sim.N1*sim.N2);

dFODC_dRho = ...
    (dFODCNotNormalized_dRho    * C - FODCNotNormalized.*dC_dRho)/C^2;
dFODC_dDelta = ...
    (dFODCNotNormalized_dDelta  * C - FODCNotNormalized.*dC_dDelta)/C^2;
dFODC_dTheta = ...
    (dFODCNotNormalized_dTheta  * C - FODCNotNormalized.*dC_dTheta)/C^2;
dFODC_dEpsilon = ...
    (dFODCNotNormalized_dEpsilon* C - FODCNotNormalized.*dC_dEpsilon)/C^2;

%% Nonlinearity derivative
ds_dodcLin = (2/omega) * s_hat .* (1 - s_hat);
ds_domega =  (-odcLin/omega).*ds_dodcLin;

%% bold and MRI Filter derivatives (to be used for energy gradients)
dFMTF_dFwhm = (-4 * FMTF.*sim.mrikrsq*pi^2*fwhm)/(8*log(2));

%% energy gradients
% prior gradients
dEprior_dRho    = dRho/rho_sigma^2;
dEprior_dDelta  = dDelta/delta_sigma^2;
dEprior_dEpsilon= dEpsilon/epsilon_sigma^2;
dEprior_dNoiseij= noise/noise_sigma^2;

% noise % FFTs AND IFFTs switched relative to TeX document!!!
% this can be done because of symmetries
dF      = fft2(d);
dFMri   = zeros(sim.N1,sim.N2);
dFMri(sim.MRIIndices1, sim.MRIIndices2) ...
    = FMTF.* dF;
dE_dNoiseij = sim.upSampleFactor^2 *ifft2(FODC.*fft2(ds_dodcLin.*...
    ifft2(makeFFT2Symmetric(dFMri))))/sigmaData^2 + dEprior_dNoiseij;

% omega
ds_domegaF = fft2(ds_domega);
dE_dOmegaToSum = ...
    ifft2(FMTF.*...
    makeFFT2Symmetric(ds_domegaF(sim.MRIIndices1,sim.MRIIndices2)));
dE_dOmega = sum(d(:).*dE_dOmegaToSum(:))/sigmaData^2;

% rho
dRhoMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dRho));
dE_dRhoToSum = ifft2(FMTF.*...
    makeFFT2Symmetric(dRhoMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE_dRho = sum(d(:).*dE_dRhoToSum(:))/sigmaData^2 + dEprior_dRho;

% delta
dDeltaMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dDelta));
dE_dDeltaToSum = ifft2(FMTF.*...
    makeFFT2Symmetric(dDeltaMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE_dDelta = sum(d(:).*dE_dDeltaToSum(:))/sigmaData^2 + dEprior_dDelta;

% epsilon
dEpsilonMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dEpsilon));
dE_dEpsilonToSum = ifft2(FMTF.*...
    makeFFT2Symmetric(dEpsilonMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE_dEpsilon = sum(d(:).*dE_dEpsilonToSum(:))/sigmaData^2 + ...
    dEprior_dEpsilon;

% theta
dThetaMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dTheta));
dE_dThetaToSum = ifft2(FMTF.*...
    makeFFT2Symmetric(dThetaMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE_dTheta = sum(d(:).*dE_dThetaToSum(:))/sigmaData^2;

% fwhm
dE_dFwhmToSum = ifft2(dFMTF_dFwhm.*...
    makeFFT2Symmetric(odcF(sim.MRIIndices1,sim.MRIIndices2)));
dE_dFwhm = sum(d(:).*dE_dFwhmToSum(:))/sigmaData^2;
end
