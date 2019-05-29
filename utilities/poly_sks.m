function [out] = poly_sks(i0,projA,projB,projC,im,ang,specData,scatParam,nPad,cg,ig,gamma)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Estimates scatter given polyquant information during an iteration
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters
% ----------
% i0            -- the incident intensity (flat field).
% projA         -- projection from first linear fit.
% projB         -- projection from second linear fit.
% projC         -- offset projection from second linear fit.
% im            -- the current estimate.
% ang           -- an array of projection angles.
% specData      -- spectrum data array.
% scatParam     -- scatter kernel parameters.
% nPad          -- padding for FFT filtering.
% cg            -- the system geometry (from Fessler's toolbox)
% ig            -- the image geometery (from Fessler's toolbox).
% gamma         -- edge compensation fudge factors:
%      gamma(1) -> edge attenuation scaling factor.
%      gamma(2) -> sets upper limit for attenuating effect.
%      gamma==0 -> no edge compensation.
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
arr1 = (-(cg.ns-1)/2-nPad:(cg.ns-1)/2+nPad)*cg.ds; %
if size(i0,2) == size(i0,1)
    arr2 = arr1;
else
    arr2 = (-(cg.nt-1)/2-nPad:(cg.nt-1)/2+nPad)*cg.dt;
end
[us,vs] = ndgrid(arr1,arr2);
uu = ((-(cg.ns-1)/2:(cg.ns-1)/2)+cg.offset_s)*cg.ds;
eDen = projA+projB;

%% Calculate the magnification factor
projScat = zeros(size(i0));
centPo = round(size(im)/2);
RP = regionprops(double(im(:,:,centPo(3))>0.1),'Centroid','MajorAxisLength','MinorAxisLength','Orientation');

tranF = (RP.Centroid-centPo([2,1]))*ig.dx;
trans = [tranF,0];
magFactor = zeros(size(i0,3),1);
tmpFac = zeros(size(i0));
for k = 1:size(i0,3)
    rotAng = ang(k);
    out = ell_centroid(rotAng,RP.MinorAxisLength*ig.dx/2,RP.MajorAxisLength*ig.dx/2,...
                  tranF(1),tranF(2),uu,true,cg);
    tmpFac(:,:,k) = repmat((cg.dod-out')./cg.dod,1,size(i0,2));
    rotAng = ang(k);
    rotMat = [cosd(rotAng),-sind(rotAng),0;
             sind(rotAng),cosd(rotAng),0;
             0,0,1];
    rotTrans = rotMat*trans';
    
    magFactor(k) = (cg.dod-rotTrans(1))/cg.dod;
    tmpFac(:,:,k) = repmat(magFactor(k),size(i0,1),size(i0,2)).^2;
    
end

%% Calculate the convolutional scatters
scat1 = zeros(size(i0,1)+2*nPad,size(i0,2)+2*nPad,size(i0,3));
g = scat1;
broad = zeros(size(i0));
specProb = specData.spectrum(:).*specData.response(:)./sum(specData.spectrum(:).*specData.response(:));

for k = 1:length(specData.energy)
    % Estimate the attenuation from the polyquant projections
    atten = specData.knee(1,1,k)*projA+specData.knee(1,2,k)*projB+specData.knee(2,2,k)*projC;
    % Narrow scatter field calculation
    A = scatParam.fA1(k,1)*ones(size(i0));...
        C = scatParam.C1(k);
    for i = 1:size(i0,3)
        A(:,:,i) = A(:,:,i)./tmpFac(:,:,i);
        tmpC = magFactor(i)*C;
        g(:,:,i) = exp(-(us.^2+vs.^2)/(tmpC^2));
    end
    
    fg = fft2(g);
    
    projFor = specProb(k)*A.*i0.*exp(-atten).*(eDen);
    scat1 = scat1+fft2(padarray(projFor,[nPad,nPad,0])).*fg;

    % Broad scatter field calculation
    A = scatParam.fA2(k,1)*ones(size(i0)); x2 = scatParam.fA2(k,2); x3 = scatParam.fA2(k,3);
    for i = 1:size(i0,3)
        A(:,:,i) = A(:,:,i)./tmpFac(:,:,i);
    end
    broadTmp = specProb(k)*A.*i0.*exp(-atten*x2).*(eDen.^x3);

    broad = broad+broadTmp;
end

if (gamma)
    broad = edge_factor(scatParam.eFac.*gamma(1),eDen,broad,abs(gamma(2)),cg);
end

gBroad = zeros(size(scat1));
for i = 1:size(i0,3)
    tmpC = sqrt(magFactor(i))*scatParam.C2(1);
    gBroad(:,:,i) = exp(-(us.^2+vs.^2)/(tmpC^2));
end

fg = fft2(gBroad);
scat1 = scat1+fft2(padarray(broad,[nPad,nPad,0])).*fg;

for j = 1:size(projA,3)
        projScat(:,:,j) = real(unpad(ifftshift(ifft2(scat1(:,:,j))),nPad));
end

projScat(projScat<0) = 0;  % ensure the scatter is non-negative

out = projScat;
end

function out = unpad(in,nPad)
sz1 = size(in,1); sz2 = size(in,2);
out = in(1+nPad:sz1-nPad,1+nPad:sz2-nPad,:);
end

function out = edge_factor(factor,eDen,broadIn,sig,cg)
%% The edge compensation factor
eDen = imgaussfilt(eDen,10);
[gradX,gradY] = gradient_op3d(eDen);

shiftFacX = eDen.*gradX.*factor./cg.ds;
shiftFacY = eDen.*gradY.*factor./cg.dt;
shiftFacX(abs(shiftFacX)>sig) = sig;
shiftFacY(abs(shiftFacX)>sig) = sig;

out = broadIn.*exp(-(shiftFacX.^2)./(35.^2)-(shiftFacY.^2)./(35.^2));
end
