function [odc,mri1,mri2,E,...
    dE_dNoiseij,dE_dRho,dE_dDelta,dE_dEpsilon,...
    dE_dTheta,dE_dOmega,dE_dFwhm1,dE_dFwhm2] = ...
    computeODCModelSimultaneous(sim,...
    noise,rho,delta,epsilon,theta,omega,beta1,beta2,fwhm1,fwhm2,...
    data1,data2,sigmaData1,sigmaData2,roi1,roi2,...
    filterCutoffs1,filterCutoffs2,limitVector)

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
sigmaPSF1 = fwhm1/(2*sqrt(2*log(2)));
sigmaPSF2 = fwhm2/(2*sqrt(2*log(2)));
FMTF1 = beta1 * exp(-2*pi^2*sigmaPSF1^2*sim.mrikrsq)/sim.upSampleFactor^2;
FMTF2 = beta2 * exp(-2*pi^2*sigmaPSF2^2*sim.mrikrsq)/sim.upSampleFactor^2;

if exist('filterCutoffs1','var') && ~isempty(filterCutoffs1)
    FMTF1((sim.mrikrsq<max(filterCutoffs1(1),0)^2)|...
        (sim.mrikrsq>filterCutoffs1(2)^2)) = 0;
end
if exist('filterCutoffs2','var') && ~isempty(filterCutoffs2)
    FMTF2((sim.mrikrsq<max(filterCutoffs2(1),0)^2)|...
        (sim.mrikrsq>filterCutoffs2(2)^2)) = 0;
end

%% compute BOLD and MRI part
odcF   = fft2(odc);
mri1    = ifft2(FMTF1.*...
    makeFFT2Symmetric(odcF(sim.MRIIndices1,sim.MRIIndices2)));
mri2    = ifft2(FMTF2.*...
    makeFFT2Symmetric(odcF(sim.MRIIndices1,sim.MRIIndices2)));

if nargout == 3
    return;
end

%% compute E
inROI1 = (roi1(:) == 1);
inROI2 = (roi2(:) == 1);
dRho    = rho - rho_mean;
dDelta  = delta - delta_mean;
dEpsilon= epsilon - epsilon_mean;
d1       = (mri1 - data1).*roi1;
d2       = (mri2 - data2).*roi2;
if fwhm1<0 || fwhm2<0 || omega<0 || delta<0 || epsilon<0 || rho<0 || ...
        rho<rhoLimit(1) || rho>rhoLimit(2) || ...
        omega<omegaLimit(1) || omega>omegaLimit(2);
    E = Inf;
else
    rhoE    = dRho^2/(2*rho_sigma^2);
    deltaE  = dDelta^2/(2*delta_sigma^2);
    epsilonE= dEpsilon^2/(2*epsilon_sigma^2);
    noiseE  = sum(noise(:).^2)./(2*noise_sigma^2);
    dataE1  = sum(d1(inROI1).^2)./(2*sigmaData1^2);
    dataE2  = sum(d2(inROI2).^2)./(2*sigmaData2^2);
    E = dataE1 + dataE2 + noiseE + rhoE + deltaE + epsilonE;
end

if nargout == 4
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
dFMTF1_dFwhm1 = (-4 * FMTF1.*sim.mrikrsq*pi^2*fwhm1)/(8*log(2));
dFMTF2_dFwhm2 = (-4 * FMTF2.*sim.mrikrsq*pi^2*fwhm2)/(8*log(2));

%% energy gradients
% prior gradients
dEprior_dRho    = dRho/rho_sigma^2;
dEprior_dDelta  = dDelta/delta_sigma^2;
dEprior_dEpsilon= dEpsilon/epsilon_sigma^2;
dEprior_dNoiseij= noise/noise_sigma^2;

% noise % FFTs AND IFFTs switched relative to TeX document!!!
% this can be done because of symmetries
dF1      = fft2(d1);
dFMri1   = zeros(sim.N1,sim.N2);
dFMri1(sim.MRIIndices1, sim.MRIIndices2) ...
    = FMTF1.* dF1;
dE1_dNoiseij = sim.upSampleFactor^2 *ifft2(FODC.*fft2(ds_dodcLin.*...
    ifft2(makeFFT2Symmetric(dFMri1))))/sigmaData1^2;

dF2      = fft2(d2);
dFMri2   = zeros(sim.N1,sim.N2);
dFMri2(sim.MRIIndices1, sim.MRIIndices2) ...
    = FMTF2.* dF2;
dE2_dNoiseij = sim.upSampleFactor^2 *ifft2(FODC.*fft2(ds_dodcLin.*...
    ifft2(makeFFT2Symmetric(dFMri2))))/sigmaData2^2;
dE_dNoiseij = dE1_dNoiseij + dE2_dNoiseij + dEprior_dNoiseij;

% omega
ds_domegaF = fft2(ds_domega);
dE1_dOmegaToSum = ...
    ifft2(FMTF1.*...
    makeFFT2Symmetric(ds_domegaF(sim.MRIIndices1,sim.MRIIndices2)));
dE1_dOmega = sum(d1(:).*dE1_dOmegaToSum(:))/sigmaData1^2;
dE2_dOmegaToSum = ...
    ifft2(FMTF2.*...
    makeFFT2Symmetric(ds_domegaF(sim.MRIIndices1,sim.MRIIndices2)));
dE2_dOmega = sum(d2(:).*dE2_dOmegaToSum(:))/sigmaData2^2;
dE_dOmega = dE1_dOmega + dE2_dOmega;

% rho
dRhoMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dRho));
dE1_dRhoToSum = ifft2(FMTF1.*...
    makeFFT2Symmetric(dRhoMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE1_dRho = sum(d1(:).*dE1_dRhoToSum(:))/sigmaData1^2;
dE2_dRhoToSum = ifft2(FMTF2.*...
    makeFFT2Symmetric(dRhoMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE2_dRho = sum(d2(:).*dE2_dRhoToSum(:))/sigmaData2^2;
dE_dRho = dE1_dRho + dE2_dRho + dEprior_dRho;

% delta
dDeltaMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dDelta));
dE1_dDeltaToSum = ifft2(FMTF1.*...
    makeFFT2Symmetric(dDeltaMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE1_dDelta = sum(d1(:).*dE1_dDeltaToSum(:))/sigmaData1^2;
dE2_dDeltaToSum = ifft2(FMTF2.*...
    makeFFT2Symmetric(dDeltaMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE2_dDelta = sum(d2(:).*dE2_dDeltaToSum(:))/sigmaData2^2;
dE_dDelta = dE1_dDelta + dE2_dDelta + dEprior_dDelta;

% epsilon
dEpsilonMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dEpsilon));
dE1_dEpsilonToSum = ifft2(FMTF1.*...
    makeFFT2Symmetric(dEpsilonMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE1_dEpsilon = sum(d1(:).*dE1_dEpsilonToSum(:))/sigmaData1^2 + ...
    dEprior_dEpsilon;
dE2_dEpsilonToSum = ifft2(FMTF2.*...
    makeFFT2Symmetric(dEpsilonMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE2_dEpsilon = sum(d2(:).*dE2_dEpsilonToSum(:))/sigmaData2^2 + ...
    dEprior_dEpsilon;
dE_dEpsilon = dE1_dEpsilon + dE2_dEpsilon;

% theta
dThetaMTFTerm = fft2(ds_dodcLin.*ifft2(noiseF.*dFODC_dTheta));
dE1_dThetaToSum = ifft2(FMTF1.*...
    makeFFT2Symmetric(dThetaMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE1_dTheta = sum(d1(:).*dE1_dThetaToSum(:))/sigmaData1^2;
dE2_dThetaToSum = ifft2(FMTF2.*...
    makeFFT2Symmetric(dThetaMTFTerm(sim.MRIIndices1,sim.MRIIndices2)));
dE2_dTheta = sum(d2(:).*dE2_dThetaToSum(:))/sigmaData2^2;
dE_dTheta = dE1_dTheta + dE2_dTheta;

% fwhm
dE_dFwhm1ToSum = ifft2(dFMTF1_dFwhm1.*...
    makeFFT2Symmetric(odcF(sim.MRIIndices1,sim.MRIIndices2)));
dE_dFwhm2ToSum = ifft2(dFMTF2_dFwhm2.*...
    makeFFT2Symmetric(odcF(sim.MRIIndices1,sim.MRIIndices2)));
dE_dFwhm1 = sum(d1(:).*dE_dFwhm1ToSum(:))/sigmaData1^2;
dE_dFwhm2 = sum(d2(:).*dE_dFwhm2ToSum(:))/sigmaData2^2;
end
