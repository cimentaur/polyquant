function [eden,mden] = mat_to_den(attenuationDb,matIm)
comp = attenuationDb.comp(1:53,:);
dens = attenuationDb.density(1:53);
Z = [1,6,7,8,11,12,15,16,17,19,20,26,53];
A = [1.01,12.01,14.01,16.00,22.99,24.305,30.97,32.066,35.45,39.098,40.08,55.845,126.90];
eDen = comp*(Z'./A');
waterE = 0.1119*1/1.01+0.8881*8/16;
relDen = dens.*eDen./waterE;
eden = matIm;
mden = matIm;
for k = 1:53
   eden(matIm==k) = relDen(k);
   mden(matIm==k) = dens(k);
end
mden(matIm==54) = 4.506;
eden(matIm==54) = 3.7326;
end