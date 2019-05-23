function out = polyquant(mode,specData,y,I0,Af,xTrue)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Performs direct quantitative reconstruction from polyergetic data.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters
% ----------
% mode          -- structure describing the dimensions of the geometry.
% specData      -- 3D volume of specimen to be projected.
% y             -- the CT measurements
% I0            -- the incident flux profile
% (xTrue)       -- ground truth image (optional)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Created:      07/03/2018
% Last edit:    26/04/2019
% Jonathan Hugh Mason
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% References: (please cite if making use of this code or its methods) 
% Jonathan H Mason et al 2017 Phys. Med. Biol. 62 8739
% Jonathan H Mason et al 2018 Phys. Med. Biol. 63 225001
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Initialisation
mode = initialise_mode(mode);
if ismatrix(y)
    x0 = ones(size(Af.arg.mask));
else
    x0 = ones(Af.arg.ig.dim);
end
Ab = Gblock(Af,mode.nSplit);

A = @(x,ind) Ab{ind}*x;
At = @(p,ind) Ab{ind}'*p;

if mode.offset
    w = @(z) offset_weight(z,mode.cg);
else
    w = @(z) z;
end
%const = y-y.*log(y+eps);
%const = sum(const(:));
if isfield(specData,'response')
    specData.spectrum = specData.spectrum.*specData.response;
end


if ~isfield(mode,'L')
    mode.L = lipscitz_estimate(specData,I0,mode.scat,y,Ab*x0,Af);
end
%
alpha = mode.nSplit*mode.tau/mode.L;  % the step-size

x1 = x0;
timeTot = tic;

if mode.nest 
   t = 1;
end

%out.f(1) = f(fShow,x1);
out.res(1) = rms(x1(:)-xTrue(:));
%if mode.verbose>0
%    fprintf('Intial objective function value = %d\n',out.f(1)-const);
%end
grAx = @(x1,is,ys,ind,subSet) polyquant_grad(specData,A,At,is,x1,ys,ind,mode.scatFun,subSet,w);

%% The main iterative loop
for k = 1:mode.maxIter
    ind = mod(k,mode.nSplit)+1;
    if mode.bitRev
       ind = bit_rev(ind-1,mode.nSplit)+1; 
    end
    subSet = ind:mode.nSplit:size(y,ndims(x0));

    if ndims(x0) == 3
        is = I0(:,:,subSet);
        ys = y(:,:,subSet);
    else
        is = I0(:,subSet);
        ys = y(:,subSet);
    end
    
    gradAx = grAx(x1,is,ys,ind,subSet);
    out.scat(:,:,subSet) = gradAx.s;
    xNew = mode.proxFun(x1-alpha*gradAx.grad,alpha);

    t1 = 0.5*(1+sqrt(1+4*t^2));
    x1 = xNew+(t-1)/t1*(xNew-x0);
    x0 = xNew;
    t = t1;

    out.res(k+1) = rms(x1(:)-xTrue(:));
    if mode.verbose > 0
      fprintf('\r res = %d        ',out.res(k+1));
    end
    if mode.verbose == 2
      fprintf(' Iter = %d; subset = %d',k,ind);
        if ndims(x1) == 3
            tmp = flipud(x1);%imrotate(x1,-90);
            subplot(1,3,1),imshow(tmp(:,:,20),[0.8,1.2]);
            subplot(1,3,2),imshow(tmp(:,:,30),[0.8,1.2]);
            subplot(1,3,3),imshow(tmp(:,:,40),[0.8,1.2]);
        else
            imshow(x1,[0,2]);
        end
        drawnow;
    end

end
time = toc(timeTot);
out.time = time;
%out.fTime = TF;
out.rec = xNew;
if mode.verbose > 0
    fprintf('finito in %d s\n',time);
end

end

function strOut = polyquant_grad(specData,A,At,I0,rho,y,ind,scatFun,subSet,w)
% This function calculates the gradient, objective function unless using
% OS, and the scatter if calculated on the fly.
projSet = cell(length(specData.hinge)-1,2);
mask = cell(length(specData.hinge)-1,1);
projSet{1,2} = 0;
for k = 1:length(specData.hinge)-1
    mask{k} = double(rho > specData.hinge(k) & rho < specData.hinge(k+1));
    projSet{k,1} = A(mask{k}.*rho,ind);
    if k>1
        projSet{k,2} = A(mask{k},ind);
    end
end
specProb = specData.spectrum./sum(specData.spectrum(:));
    
mainFac = zeros(size(y));
hingeFac = cell(length(specData.hinge)-1);
for k = 1:length(specData.hinge)-1
    hingeFac{k} = zeros(size(y));
end

s = scatFun(I0,projSet{1,1},projSet{2,1},projSet{2,2},rho,subSet,specData.knee);
for k = 1:length(specData.spectrum)
    linSum = zeros(size(y));
    for l = 1:length(specData.hinge)-1
        linSum = linSum+specData.knee(1,l,k)*projSet{l,1}...
                       +specData.knee(2,l,k)*projSet{l,2};
    end
    tmp = specProb(k).*exp(-linSum);
    mainFac = mainFac+tmp;
    for l = 1:length(specData.hinge)-1
        hingeFac{l} = hingeFac{l}+tmp*specData.knee(1,l,k);
    end
end
mainFac = I0.*mainFac;

deriFac = w(y./(mainFac+s)-1);

out = zeros(size(rho));
for l = 1:length(specData.hinge)-1
    out = out+mask{l}.*At(I0.*hingeFac{l}.*deriFac,ind);
end

strOut.grad = out;
strOut.s = s;
end

function out = lipscitz_estimate(specData,I0,s,y,flat,At)
% A crude but reasonably acceptable estimate of the Lipschitz constant
specProb = specData.spectrum./sum(specData.spectrum(:));
tmpA = 0;
for k = 1:length(specData.spectrum)
    tmpA = tmpA+specProb(k)*specData.knee(1,1,k).^2;
end
fac = I0.*(1-y.*s./((I0+s).^2));
p2A = At'*(flat.*tmpA.*fac);
out = max(p2A(:));
end

function out = prox_nz(in,up)
% Simple proximal function to enforce box constraints
if nargin > 1
    in(in>up) = up;
end
in(in<0) = 0;
out = in;
end

function mode = initialise_mode(mode)
% Make sure everything is in order
if ~isfield(mode,'nest'),       mode.nest = true; end
if ~isfield(mode,'maxIter'),    mode.maxIter = 100; end
if ~isfield(mode,'bitRev'),     mode.bitRev = true; end
if ~isfield(mode,'offset'),     mode.offset = false; end
if ~isfield(mode,'verbose'),    mode.verbose = 1; end
if ~isfield(mode,'tau'),        mode.tau = 1.99; end
if ~isfield(mode,'nSplit'),     mode.nSplit = 1; end
if ~isfield(mode,'flip'),       mode.flip = false; end
if ~isfield(mode,'regFun'),     mode.regFun = @(z) 0; end
if ~isfield(mode,'proxFun'),    mode.proxFun = @(z,t) prox_nz(z); end
if ~isfield(mode,'scatFun')
    mode.scat = 0;
    mode.scatFun = @(z,~,~,~,~,~,~) 0; 
elseif ~isa(mode.scatFun,'function_handle')
    mode.scat = mode.scatFun;
    mode.scatFun = @(z,~,~,~,~,subSet,~) mode.scat(:,:,subSet);
else
    mode.scat = 0;
end
end

function out = offset_weight(proj,cg)
    out = proj;
    us = ((cg.ns/2-0.5):-1:(-cg.ns/2+0.5))*cg.ds - cg.offset_s*cg.ds;
    overlap = max(us);
    overLoc = sum(abs(us)<=overlap);
    replaceLoc = 1:overLoc;
    denom = 2*atan(overlap/cg.dsd);
    num = pi*atan(us(replaceLoc)/cg.dsd);
    %weightArray = 1-cos(linspace(0,pi/2,overLoc)).^2;
    weightArray = 1-0.5*(sin(num./denom)+1);
    weightMat = repmat(weightArray',1,cg.nt);
    replaceLoc = 1:size(weightMat,1);

    for k = 1:size(proj,3)
        out(end-replaceLoc+1,:,k) = proj(end-replaceLoc+1,:,k).*weightMat;
    end
end