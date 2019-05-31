%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Fanbeam Polyquant demo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Description
% ----------
% This script demonstrates direct quantitative reconstruction into relative
% electron density from polyenergetic X-ray CT measurements. There are 6
% samples, including various anatomical sites from the ICRP 110 female
% computational phantom, and one ('implant') with added double titanium
% hips.
% The measurements were generated using Gate software and consist of:
% proj -- the total measured photons after detector response function.
% scat -- the photons from scatter (very small and ignored here).
% i0   -- the incident source flux.
% The structure in fandata.specData contains the spectral information of
% the source (only a subsampling of source used in simulation):
% specData.energy   -- the energies (MeV) in the subsampled spectrum.
% specData.spectrum -- the subsampled source spectrum.
% specData.response -- the detector response function.
% specData.hinge    -- the location of the piecewise linear fit
%                      transitions, for 3 linear sections.
% specData.knee     -- contains the equations for the piecewise linear fits
%                      between relative electron density and each energy in
%                      specData.energy. This was fitted against the
%                      biological materials in the ICRP 89 and for titanium
%                      (density = 4.506 g/cm3).
% Things to try
% ------------
% o Test the various samples by changing 'sample' string.
% o Test the influence on regularisation constant 'lambda'.
% o Try removing the 'mode.proxFun' line, to give ML reconstruction.
% o Adjust 'mode.maxIter', 'mode.tau' and 'mode.nSplit' to change
%   convergence properties.
% o Set 'mode.numLinFit = 1' to see beam-hardening artefacts, which is
%   equivalent to a linearised model (monoenergetic).
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
load data/attenuationDb
load data/fandata
addpath(genpath('.'));

%% Setup the geometry 
ig = image_geom('nx',137,'ny',299,'dx',0.1775,'dy',0.1775);
sg = sino_geom('fan','dfs',inf,'dsd',120.3,'dso',64.5,'units','cm',...
               'ns',512,'ds',85/512,'strip_width','ds','na',360,'down',1);
A = Gtomo2_dscmex(sg, ig);

%% Polyquant setup 
sample = 'chest';  % 'brain', 'head', 'abdomen', 'pelvis', 'implant'
mode = [];
mode.useConst = true;  % just to offset objective function to better range
mode.contrast = [0.5,1.4];  % display contrast
mode.verbose = 2;  % 0 = no output; 1 = text output; 2 = text+image output
mode.tau = 5;  % conservative choice
mode.nSplit = 24;
mode.maxIter = 500;  % more iterations recommended for implant
lambda = 0.5;  % can be optimised for better results
mode.numLinFit = 2;  % can be set to 2 for tissue, but use 3 for implant
mode.proxFun = @(z,t) prox_tv_nn(z,t*lambda);
mode.regFun = @(z) norm_tv(z);
xTrue = mat_to_den(attenuationDb,single(fandata.(sample).mat));

%% Perform the reconstruction
out = polyquant(mode,fandata.specData,fandata.(sample).proj,fandata.i0,A,xTrue);
fprintf('Reconstructed with PSNR = %.2f dB\n',20*log10(max(xTrue(:))./out.rmse(end)));