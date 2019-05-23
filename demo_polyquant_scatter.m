pelvis.cg = ct_geom('fan', 'ns', 256, 'nt', 128, 'na', 160, ...
        'orbit_start',-90, 'orbit',-360,... 201.528 ,...+1.125, ...
		'ds', 0.1552, 'dt', 0.2328, ...
		'offset_s', 16/0.1552, ... % quarter detector
		'offset_t', 0.0, ...
		'dsd', 150,'dod', 50, 'dfs', inf);
pelvis.ig = image_geom('nx', 299, 'ny', 137, 'nz', 60, 'dx', 0.1775,'dy',0.1775,'dz',0.484);%320);%,'mask',bshadow==1);
pelvis.A = Gcone(pelvis.cg, pelvis.ig, 'type', 'sf2', 'class', 'Fatrix');
%%
mode = [];
scatParam.ang = ((-90:2.25:270-2.25));
mode.verbose = 2;
mode.tau = 2;
mode.offset = true;
mode.cg = pelvis.cg;
mode.nSplit = 32;
mode.maxIter = 500;
mode.scatFun = @(I0,projA,projB,projC,rho,subSet,knee) poly_sks(1.5*I0,projA,projB,projC,rho,scatParam.ang(subSet),...
                   specData,scatFacFix,32,pelvis.cg,false,false,[0.3,15],1);
lambda = 1;
param.up = 1.8;
mode.proxFun = @(z,t) prox_tv3d_nn(z,t*lambda,param);
xTrue = pelvis.eden;
out = polyquant(mode,specData,flipud(pelvis.fac*pelvis.proj),flipud(i0),pelvis.A,xTrue);