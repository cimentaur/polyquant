%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cone-beam CT Polyquant with scatter model demo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Description
% ----------
% This script demonstrates polyquant reconstruction for cone-beam CT, with
% high levels of scatter. The data was generated using Gate without a
% collimator. Please refer to 'demo_polyquant_fanbeam.m' for a description
% of the variables, which are the same in this case, only with a different
% spectrum (pelvis.specData.spectrum).
% In this demo, the scatter esimation function 'poly_sks' is called at each
% iteration.
% This is data is very extreme in terms of scatter and photon starvation,
% so is challenging to generate accurate results.
%
% Things to try
% ------------
% o See the influence of the edge compensation factor by replacing the
%   '[0.3,15]' parameters with '0', which will disable it.
% o Compare the estimated scatter in 'out.scat' against the true scatter
%   from 'pelvis.scat'.
% o Substitute the scatter estimate (mode.scatFun) for the true scatter
%   (pelvis.scat), to see the ultimate scatter estimate's performance.
% o Remove the 'mode.scatFun' line to see the scatter artefacts.
% o Adjust regularisation and convergence parameters an compare to fanbeam.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Created:      26/04/2019
% Last edit:    29/05/2019
% Jonathan Hugh Mason
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% References: (please cite if making use of this code or its methods) 
% Jonathan H Mason et al 2017 Phys. Med. Biol. 62 8739
% Jonathan H Mason et al 2018 Phys. Med. Biol. 63 225001
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Load data
load data/scat_param
load data/data_cbct_pelvis
addpath(genpath('.'));
i0 = repmat(pelvis.i0,1,1,160);

%% Setup the geometry 
cg = ct_geom('fan', 'ns', 256, 'nt', 128, 'na', 160, ...
        'orbit_start',-90, 'orbit',-360,... 201.528 ,...+1.125, ...
		'ds', 0.1552, 'dt', 0.2328, ...
		'offset_s', 16/0.1552, ... % quarter detector
		'offset_t', 0.0, ...
		'dsd', 150,'dod', 50, 'dfs', inf);
ig = image_geom('nx', 299, 'ny', 137, 'nz', 60, 'dx', 0.1775,'dy',0.1775,'dz',0.484);
A = Gcone(cg, ig, 'type', 'sf2', 'class', 'Fatrix');
%% Polyquant setup and reconstruction
mode = [];
mode.verbose = 2;
mode.tau = 5;  % a little more aggressive
mode.offset = true;
mode.cg = cg;
mode.nSplit = 32;
mode.maxIter = 500;
mode.numLinFit = 2;  % since there's no metal implants
lambda = 0.1;  % can be optimised for better results
mode.proxFun = @(z,t) prox_tv3d_nn(z,t*lambda);
angArray = ((-90:2.25:270-2.25));
%% Scatter estimation function
% The factor 1.5 on the incident intensity is to compensate for the 
% detector not having square elements.
% The edge factor 0.3 can be increased to around 0.5 for full-fan scanning.
mode.scatFun = @(i0,projA,projB,projC,rho,subSet,knee) ...
                 poly_sks(1.5*i0,projA,projB,projC,rho,angArray(subSet),...
                 pelvis.specData,scatParam,32,cg,ig,[0.3,15]);
%% Perfrm Polyquant reconstrution
out = polyquant(mode,pelvis.specData,pelvis.proj,i0,A,pelvis.eden);