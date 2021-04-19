function sim = setupsim(nVoxels1,nVoxels2,vSize,upSampleFactor)
sim.upSampleFactor = upSampleFactor;
sim.nVoxels1 = nVoxels1;
sim.nVoxels2 = nVoxels2;
sim.vSize = vSize;
sim.N1 = nVoxels1 * upSampleFactor;
sim.L1 = nVoxels1 * vSize;
sim.N2 = nVoxels2 * upSampleFactor;
sim.L2 = nVoxels2 * vSize;
sim.dx = vSize/upSampleFactor;
sim.dk1 = 1/sim.L1;
sim.dk2 = 1/sim.L2;
sim.x1Range = ifftshift((-ceil((sim.N1-1)/2):floor((sim.N1-1)/2))*sim.dx);
sim.x2Range = ifftshift((-ceil((sim.N2-1)/2):floor((sim.N2-1)/2))*sim.dx);
sim.k1Range = ifftshift((-ceil((sim.N1-1)/2):floor((sim.N1-1)/2))*sim.dk1);
sim.k2Range = ifftshift((-ceil((sim.N2-1)/2):floor((sim.N2-1)/2))*sim.dk2);
[sim.k1,sim.k2] = ndgrid(sim.k1Range,sim.k2Range);
sim.krsq = sim.k1.^2 + sim.k2.^2;
sim.kr   = sqrt(sim.krsq);
sim.kphi = atan2(sim.k2,sim.k1);
if mod(nVoxels1,2) || mod(nVoxels2,2)
    error('Even number of voxels per dimension required!');
else
    NMRIhalf1 = nVoxels1/2;
    NMRIhalf2 = nVoxels2/2;
end
sim.MRIIndices1 = [1:NMRIhalf1 sim.N1-NMRIhalf1+1:sim.N1];
sim.MRIIndices2 = [1:NMRIhalf2 sim.N2-NMRIhalf2+1:sim.N2];
sim.mrikrsq = sim.krsq(sim.MRIIndices1,sim.MRIIndices2);
end
