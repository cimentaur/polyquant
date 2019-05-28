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
mode.contrast = [0.5,1.4];  % display contrast
mode.verbose = 2;
mode.tau = 2;  % conservative choice
mode.nSplit = 24;
mode.maxIter = 500;  % more iterations recommended for implant
lambda = 0.5;  % can be optimised for better results
mode.numLinFit = 2;  % can be set to 2 for tissue, but use 3 for implant
mode.proxFun = @(z,t) prox_tv_nn(z,t*lambda);
xTrue = mat_to_den(attenuationDb,single(fandata.(sample).mat));

%% Perform the reconstruction
out = polyquant(mode,fandata.specData,fandata.(sample).proj,fandata.i0,A,xTrue);