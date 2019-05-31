%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cone-beam CT Polyquant demo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Description
% ----------
% This script demonstrates polyquant reconstruction for cone-beam CT, with
% high levels of scatter. The data was generated using Gate without a
% collimator. Please refer to 'demo_polyquant_fanbeam.m' for a description
% of the variables, which are the same in this case, only with a different
% spectrum (head.specData.spectrum).
% In this demo, we are given a pre-calculated estimate of the scatter
% (scatEst) generated from an fASKS-like approach. 
%
% Things to try
% ------------
% o Substitute the scatter estimate (mode.scatFun) for the true scatter
%   (head.scat), to see the ultimate scatter estimate's performance.
% o Remove the 'mode.scatFun' line to see the scatter artefacts.
% o Adjust regularisation and convergence parameters an compare to fanbeam.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Created:      26/04/2019
% Last edit:    31/05/2019
% Jonathan Hugh Mason
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% References: (please cite if making use of this code or its methods) 
% Jonathan H Mason et al 2017 Phys. Med. Biol. 62 8739
% Jonathan H Mason et al 2018 Phys. Med. Biol. 63 225001
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Load data
load data/scatter_est_head
load data/data_cbct_head
addpath(genpath('.'));
i0 = repmat(head.i0,1,1,160);

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
mode.useConst = true;
mode.verbose = 2;
mode.tau = 2;
mode.nSplit = 32;
mode.maxIter = 300;
mode.scatFun = scatEst;
lambda = 2;  % can be optimised for better results
mode.proxFun = @(z,t) prox_tv3d_nn(z,t*lambda);
mode.regFun = @(z) norm_tv3d(z);
out = polyquant(mode,head.specData,head.proj,i0,A,head.eden);
fprintf('Reconstructed with PSNR = %.2f dB\n',20*log10(max(head.eden(:))./out.rmse(end)));