%% Setup the geometry 
ig = image_geom('nx',137,'ny',299,'dx',0.1775,'dy',0.1775);
sg = sino_geom('fan','dfs',inf,'dsd',120.3,'dso',64.5,'units','cm',...
               'ns',512,'ds',85/512,'strip_width','ds','na',360,'down',1);
A = Gtomo2_dscmex(sg, ig);

%%
[specData.knee,specData.fac] = pw_knee_fit(specData,attenuationDb,false,tiMa);
specData.hinge = [0;specData.fac;inf];
specData.kneeOld = specData.knee;
specData.knee = zeros(2,3,21);
specData.knee(1,1,:) = permute(specData.kneeOld(1,:),[3,1,2]);
specData.knee(1,2,:) = permute(specData.kneeOld(2,:),[3,1,2]);
specData.knee(2,2,:) = permute(specData.kneeOld(3,:),[3,1,2]);
specData.knee(1,3,:) = permute(specData.kneeOld(4,:),[3,1,2]);
specData.knee(2,3,:) = permute(specData.kneeOld(5,:),[3,1,2]);
%%
mode = [];
mode.verbose = 2;
mode.tau = 1e5;
mode.nSplit = 16;
mode.maxIter = 5e2;
lambda = 0.5;
mode.proxFun = @(z,t) prox_tv(z,t*lambda);
out = polyquant(mode,specData,fac*mc.chest.proj,i0,A,im.eden.chest);