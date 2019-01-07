function out = polyquant(mode,specData,ADB,y,scatFun,I0,Af,spec,lambda)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generates a forward cone-beam projection of a given specimen using 3D 
% interpolation. Can operate in single or double point precision.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters
% ----------
% geometry      -- structure describing the dimensions of the geometry.
% specimen      -- 3D volume of specimen to be projected.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Created:      07/03/2018
% Last edit:    07/01/2019
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Jonathan Hugh Mason
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tiMa = specData.tiMa;
opts.maxIters = mode.maxIters;
opts.verbose = mode.verbose;
Ab = Gblock(Af,mode.nSplit);
if mode.flip
    A = @(x,ind) (flipud(Ab{ind}*x));
    At = @(p,ind) (Ab{ind}'*flipud(p));
else
    A = @(x,ind) Ab{ind}*x;
    At = @(p,ind) Ab{ind}'*p;
end
const = y-y.*log(y+eps);
const = sum(const(:));
specData.spectrum = specData.spectrum.*specData.response;
[specData.knee,specData.fac] = pw_knee_fit(specData,ADB,0,tiMa);
x0 = spec;%ones(size(spec));
param.up = 8;
scat = zeros(size(y));
param.weights = ones(size(x0));
param.maxit = 5;
if strcmp(mode.reg,'tv')
    if ndims(x0) == 3
        f = @(z,x) z+lambda*norm_tv3d(x);
        prox = @(z,t,p) prox_tv3d(z,lambda,param);
    else
        f = @(z,x) z+double(lambda*norm_tv(x));
        prox = @(z,t,p) prox_tv(z,lambda,param);
    end
elseif strcmp(mode.reg,'huber')
    R = Reg1(ones(size(x0)), 'type_denom', 'matlab', ...
        'pot_arg', {'huber', lambda(1)}, 'beta', lambda(2));
    f = @(z,x) z+R.penal(R,x);
    prox = @(z) prox_nz(z-mode.tau*(1-0.8)/mode.L*R.cgrad(R,z),2); 
else
    f = @(z,x) z;
    prox = @(z,t,p) prox_nz(z,param.up);
end
hess.use = false;

if isfield(mode,'L')
    L = mode.L;
else
    hess.flat = Ab*x0;
    hess = hess_p2(specData,I0,zeros(size(I0)),y,hess,Af);
    L = max(hess.p2A(:));
end
%
beta = mode.beta;
alpha = mode.nSplit*mode.tau*(1-beta)/L; % 60

x1 = x0;
timeTot = tic;
if mode.offset
    w = @(z) offset_weight(z,mode.cg,0);
else
    w = @(z) z;
end
if mode.nest 
   t = 1; v = 0; tSum = t;
end
timeF = tic;
if mode.calcF
    fShow = poly_f(specData,Af,I0,x1,scatEst,y);
else
    fShow = 0;
end
TF = toc(timeF);
out.f(1) = f(fShow,x1);
out.res(1) = rms(x1(:)-spec(:));
if mode.verbose>0
    fprintf('Intial objective function value = %d\n',out.f(1)-const);
end

grAx = @(x1,is,ys,ind,subSet) poly_quantf(specData,A,At,is,x1,ys,ind,w,scatFun,subSet);

for k = 1:opts.maxIters
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
    scat(:,:,subSet) = gradAx.s;
    xNew = prox(x1-alpha*gradAx.grad+beta*(x1-x0));
    if mode.nest
       if mode.nest == 2
            v = v+t*alpha*gradAx.grad;
            t = 0.5*(1+sqrt(1+4*t^2));
            tSum = tSum+t;
            xNew = xNew+t/tSum*(prox(x0-v)-xNew);
            x1 = xNew;
       else
            t1 = 0.5*(1+sqrt(1+4*t^2));
            x1 = xNew+(t-1)/t1*(xNew-x0);
            x0 = xNew;
            t = t1;
       end
    else
        x0 = x1;
        x1 = xNew;
    end
    
      timeF = tic;
      if mode.calcF
          fShow = poly_f(specData,Af,I0,x1,scatEst,y);
      else
          fShow = 0;%
      end
      TF = TF+toc(timeF);
      out.f(k+1) = f(fShow,x1);
      out.res(k+1) = rms(x1(:)-spec(:));
      if opts.verbose > 0
          fprintf('f = %d; res = %d; diff = %d\n',out.f(k+1)-const,out.res(k+1),...
                  out.f(max(1,k))-out.f(k+1));
      end
      if opts.verbose == 2
          fprintf('Iter = %d; subset = %d\n',k,ind);
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
      if out.f(k+1)-const<mode.tol
          out.iters = k;
          out.converge = true;
          break;
      end
      if toc(timeTot)-TF>mode.maxTime
         out.timeOut = true;
         break;
      end
end
time = toc(timeTot);
out.time = time;
out.fTime = TF;
out.rec = xNew;
out.scat = scat;
if opts.verbose > 0
    fprintf('finito in %d s\n',time);
end

end

function strOut = poly_quantf(specData,A,At,I0,rho,y,ind,w,scatFun,subSet)
%specProb = specData.tweaker./sum(specData.tweaker);
% Construct knee function
x0 = specData.fac(1);
f = double(rho>x0);
projA = A((1-f).*rho,ind); projB = A(f.*rho,ind); projC = A(f,ind);
facA = (1-f);
facB = f;
if specData.bt
    if ndims(y) == 3
        specProb = repmat(specData.tweaker(:,:,1),1,1,size(y,3));
    else
        specProb = repmat(specData.tweaker(:,1),1,size(y,2));
    end
else
    specRat = specData.spectrum./sum(specData.spectrum(:));
    specProb = specRat(1);
end
tmp1 = specProb.*exp(-(specData.knee(1,1)*projA+...
                       specData.knee(2,1)*projB+specData.knee(3,1)*projC));
tmpA = tmp1*specData.knee(1,1);
tmpB = tmp1*specData.knee(2,1);
for k = 2:length(specData.spectrum)
    if specData.bt
        if ndims(y) == 3
            specProb = repmat(specData.tweaker(:,:,k),1,1,size(y,3));
        else
            specProb = repmat(specData.tweaker(:,k),1,size(y,2));
        end
    else
        specProb = specRat(k);
    end
    tmp = specProb.*exp(-(specData.knee(1,k)*projA+...
                       specData.knee(2,k)*projB+specData.knee(3,k)*projC));
    tmp1 = tmp1+tmp;
    tmpA = tmpA+tmp*specData.knee(1,k);
    tmpB = tmpB+tmp*specData.knee(2,k);
end
tmp1 = I0.*tmp1;
tmpA = I0.*tmpA;
tmpB = I0.*tmpB;
%%%%%%%%%%%%%%%%%%%%%%% Scatter calculation Zone %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = scatFun(I0,projA,projB,projC,rho,subSet,specData.knee);
deriFac = (y./(tmp1+s)-1);
outA = w(tmpA.*deriFac);
outB = w(tmpB.*deriFac);
out = facA.*At(outA,ind)+facB.*At(outB,ind);

strOut.grad = out;
strOut.s = s;
end

function out = hess_p2(specData,I0,s,y,in,At)
out = in;
%specProb = specData.spectrum./sum(specData.spectrum);
if specData.bt
    if ndims(y) == 3
        specProb = repmat(specData.tweaker(:,:,1),1,1,size(y,3));
    else
        specProb = repmat(specData.tweaker(:,1),1,size(y,2));
    end
else
    specRat = specData.spectrum./sum(specData.spectrum(:));
    specProb = specRat(1);
end
tmpA = specProb*specData.knee(1,1).^2;
tmpB = specProb*specData.knee(2,1).^2;
for k = 2:length(specData.spectrum)
    if specData.bt
        if ndims(y) == 3
            specProb = repmat(specData.tweaker(:,:,k),1,1,size(y,3));
        else
            specProb = repmat(specData.tweaker(:,k),1,size(y,2));
        end
    else
        specProb = specRat(k);
    end
    tmpA = tmpA+specProb*specData.knee(1,k).^2;
    tmpB = tmpB+specProb*specData.knee(2,k).^2;
end
fac = I0.*(1-y.*s./((I0+s).^2));
out.p2A = At'*(in.flat.*tmpA.*fac);
out.p2B = At'*(in.flat.*tmpB.*fac);
end

function out = prox_nz(in,up)
in(in>up) = up;
in(in<0) = 0;
out = in;
end