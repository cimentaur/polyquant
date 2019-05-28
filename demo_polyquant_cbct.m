%% Load data
load data/data_cbct_head
addpath(genpath('.'));

%% Setup the geometry 
cg = ct_geom('fan', 'ns', 256, 'nt', 128, 'na', 160, ...
        'orbit_start',-90, 'orbit',-360,...
		'ds', 0.1552, 'dt', 0.2328, ...
		'offset_s', 0, ...
		'offset_t', 0.0, ...
		'dsd', 150,'dod', 50, 'dfs', inf);
ig = image_geom('nx', 99, 'ny', 137, 'nz', 70, 'dx', 0.1775,'dy',0.1775,'dz',0.484);
A = Gcone(cg, ig, 'type', 'sf2', 'class', 'Fatrix');

%% Polyquant setup and reconstruction
mode = [];
mode.verbose = 2;
mode.tau = 2;
mode.nSplit = 32;
mode.maxIter = 300;
mode.scatFun = head.scatEst;
lambda = 2;  % can be optimised for better results
mode.proxFun = @(z,t) prox_tv3d_nn(z,t*lambda);
out = polyquant(mode,head.specData,head.proj,head.i0,A,head.eden);