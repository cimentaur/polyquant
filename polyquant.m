function out = polyquant(mode,specData,y,I0,Af,xTrue)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Performs direct quantitative reconstruction from polyergetic data.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters
% ----------
% mode          -- structure containing the settings and functions:
% (all these settings have default values: see initialise_mode)
%   mode.tau          -- stepsize scaling factor (< 2 is conservative).
%   mode.maxIter      -- number of iterations.
%   mode.nest         -- use FISTA-like Nesterov acceleration.
%   mode.nSplit       -- number of ordered subset divisions (1 is full).
%   mode.verbose      -- output settings: 0 = silent; 1 = text; 2 = figure.
%   mode.contrast     -- display contrast for output live updat figure.
%   mode.regFun       -- handle to regularisation function.
%   mode.proxFun      -- handle to proximity operator for regularisation.
%   mode.scatFun      -- scatter estimation function (see poly_sks.m).
%   mode.useConst     -- offset objective function to better range.
%   mode.bitRev       -- use subset shuffling (bit-reversal ordering).
%   mode.offset       -- use Wang offset detector weighting for half-fan.
%   mode.L            -- supplying Lipschitz estimate will save time.
% specData      -- structure containing spectral information:
%   specData.energy   -- the energies (MeV) in the subsampled spectrum.
%   specData.spectrum -- the subsampled source spectrum.
%   specData.response -- the detector response function.
%   specData.hinge    -- the location of the piecewise linear fit
%                        transitions, for 3 linear sections.
%   specData.knee     -- contains the equations for the piecewise linear 
%                        fits between relative electron density and each 
%                        energy in specData.energy. This was fitted against 
%                        the biological materials in the ICRP 89 and for
%                        titanium (density = 4.506 g/cm3).
% y             -- the raw X-ray CT measurements.
% I0            -- the incident flux profile.
% Af            -- the CT system operator generated from Fessler's toolbox. 
% xTrue         -- ground truth image (can be 0 if unknown).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Created:      07/03/2018
% Last edit:    02/06/2019
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

if isfield(mode,'numLinFit')
    specData.hinge = [specData.hinge(1:mode.numLinFit);inf];
end

if isfield(specData,'response')
    specData.spectrum = specData.spectrum.*specData.response;
end


if ~isfield(mode,'L')  % estimate Lipschitz if unknown
    mode.L = lipscitz_estimate(specData,I0,mode.scat,y,Ab*x0,Af);
end

alpha = mode.nSplit*mode.tau/mode.L;  % the step-size
if mode.useConst
    const = y-y.*log(y+eps);
    const = sum(const(:));  % a constant offset for objective function
else
    const = 0;
end

x1 = x0;
timeTot = tic;

if mode.nest 
   t = 1;
end

out.rmse(1) = rms(x1(:)-xTrue(:));
if mode.verbose == 2
    if ndims(xTrue) == 3
        subplot(2,3,1),imshow(imrotate(xTrue(:,:,20),-90),mode.contrast);
        subplot(2,3,2),imshow(imrotate(xTrue(:,:,30),-90),mode.contrast),title('ground truth');
        subplot(2,3,3),imshow(imrotate(xTrue(:,:,40),-90),mode.contrast);
    else
        subplot(2,1,1),imshow(xTrue,mode.contrast),title('ground truth');
        subplot(2,1,2)
    end
    drawnow;
end
grAx = @(x1,is,ys,ind,subSet) polyquant_grad(specData,A,At,is,x1,ys,ind,mode.scatFun,subSet,w);
objFac = zeros(size(y)); out.scat = zeros(size(y));
%% The main iterative loop
if mode.verbose > 0
    fprintf('Starting Polyquant reconstruction:\n');
end
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
    if ndims(x0) == 3
        out.scat(:,:,subSet) = gradAx.s;
        objFac(:,:,subSet) = gradAx.objFac;
    else
        out.scat(:,subSet) = gradAx.s;
        objFac(:,subSet) = gradAx.objFac;
    end
    xNew = mode.proxFun(x1-alpha*gradAx.grad,alpha);
    
    if mode.nest
        t1 = 0.5*(1+sqrt(1+4*t^2));
        x1 = xNew+(t-1)/t1*(xNew-x0);
        x0 = xNew;
        t = t1;
    else
        x1 = xNew;
    end
    
    out.rmse(k+1) = rms(x1(:)-xTrue(:));
    out.obj(k+1) = sum(double(objFac(:)+out.scat(:)-y(:).*log(objFac(:)+out.scat(:)+eps)))-const+mode.regFun(x1);
    if mode.verbose > 0
      fprintf('\rIter = %i;\t RMSE = %.4e;\t obj = %.4e;\t subset = %i    ',k,out.rmse(k+1),out.obj(k+1),ind);
    end
    if mode.verbose == 2
        str = ['polyquant at iteration: ',num2str(k)];
        if ndims(x1) == 3
            subplot(2,3,4),imshow(imrotate(x1(:,:,20),-90),mode.contrast);
            subplot(2,3,5),imshow(imrotate(x1(:,:,30),-90),mode.contrast),title(str);
            subplot(2,3,6),imshow(imrotate(x1(:,:,40),-90),mode.contrast);
        else
            imshow(x1,mode.contrast),title(str);
        end
        drawnow;
    end

end
time = toc(timeTot);
out.time = time;
out.recon = xNew;
if mode.verbose > 0
    fprintf('\n Finished in %.2e seconds\n',time);
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

if length(specData.hinge)>2  % to bodge error for one linear fit
    s = scatFun(I0,projSet{1,1},projSet{2,1},projSet{2,2},rho,subSet,specData.knee);
else
    s = scatFun(I0,projSet{1,1},projSet{1,1},projSet{1,2},rho,subSet,specData.knee);
end
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
strOut.objFac = mainFac;
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
if ~isfield(mode,'contrast'),   mode.contrast = [0,2]; end
if ~isfield(mode,'useConst'),   mode.useConst = false; end
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
    % Offset weighting for half-fan case from [G. Wang, Med Phys. 2002]
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