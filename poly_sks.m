function [out] = poly_sks(i0,projA,projB,projC,im,ang,specData,scatFac,nPad,cg,collim,bt,gamma,val)
arr1 = (-(cg.ns-1)/2-nPad:(cg.ns-1)/2+nPad)*cg.ds; %
if size(i0,2) == size(i0,1)
    arr2 = arr1;
else
    arr2 = (-(cg.nt-1)/2-nPad:(cg.nt-1)/2+nPad)*cg.dt;
end
[us,vs] = ndgrid(arr1,arr2);
uu = ((-(cg.ns-1)/2:(cg.ns-1)/2)+cg.offset_s)*cg.ds;
eDen = projA+projB;

% Calculate the magnification factor
projScat = zeros(size(i0));
centPo = round(size(im)/2);
%stats = regionprops(ones(size(im)),im,'WeightedCentroid','Centroid');
RP = regionprops(double(im(:,:,centPo(3))>0.1),'Centroid','MajorAxisLength','MinorAxisLength','Orientation');

%trans = (stats.WeightedCentroid-stats.Centroid)*0.1775;
tranF = (RP.Centroid-centPo([2,1]))*0.1776;
trans = [tranF,0];
magFactor = zeros(size(i0,3),1);
tmpFac = zeros(size(i0));
for k = 1:size(i0,3)
    rotAng = ang(k);
    [out,flTmp] = ell_centroid(rotAng,RP.MinorAxisLength*0.1775/2,RP.MajorAxisLength*0.1775/2,...
                  tranF(1),tranF(2),uu,true);
    tmpFac(:,:,k) = repmat((50-out')./50,1,size(i0,2));
    rotAng = ang(k);
    rotMat = [cosd(rotAng),-sind(rotAng),0;
             sind(rotAng),cosd(rotAng),0;
             0,0,1];
    rotTrans = rotMat*trans';
    if val == 1
        magFactor(k) = (50-rotTrans(1))/50;
        tmpFac(:,:,k) = repmat(magFactor(k),size(i0,1),size(i0,2)).^2;
    elseif val == 2
        magFactor(k) = (50-out(25))/50;
    elseif val == 3
        magFactor(k) = (50-out(25))/50;
        tmpFac(:,:,k) = repmat(magFactor(k),256,128).^2;
    elseif val == 4
        magFactor(k) = (50-mean(out(~flTmp)))/50;
        %tmpFac(:,:,k) = repmat(magFactor(k),256,128).^2;
    else
        magFactor(k) = 1;
        tmpFac(:,:,k) = 1;
    end
end
if val == 2
    tmpFac = tmpFac.^2;
end
triang = 1-abs(arr2)*7/50;
triang(triang<=0.3) = 0.3;
triang = repmat(triang,length(arr1),1);
scat1 = zeros(size(i0,1)+2*nPad,size(i0,2)+2*nPad,size(i0,3));
g = scat1;
broad = zeros(size(i0));
specProb = specData.spectrum(:).*specData.response(:)./sum(specData.spectrum(:).*specData.response(:));

for k = 1:length(specData.energy)
    atten = specData.knee(1,1,k)*projA+specData.knee(1,2,k)*projB+specData.knee(2,2,k)*projC;
    % Narrow scatter field calculation
    A = scatFac.fA1(k,1)*ones(size(i0));...
        C = scatFac.C1(k);
    for i = 1:size(i0,3)
        %A(:,:,i) = A(:,:,i)/magFactor(i)^2;
        A(:,:,i) = A(:,:,i)./tmpFac(:,:,i);
        tmpC = magFactor(i)*C;
        g(:,:,i) = exp(-(us.^2+vs.^2)/(tmpC^2));
    end
    if collim
        fg = fft2(g.*triang);
    else
        fg = fft2(g);
    end
    if bt
        specProb = specData.tweaker(:,:,k);
        projFor = specProb.*A.*i0.*exp(-atten).*(eDen);
    else
        projFor = specProb(k)*A.*i0.*exp(-atten).*(eDen);
    end
    scat1 = scat1+fft2(padarray(projFor,[nPad,nPad,0])).*fg;

    %% Broad scatter field calculation
    A = scatFac.fA2(k,1)*ones(size(i0)); x2 = scatFac.fA2(k,2); x3 = scatFac.fA2(k,3);
    for i = 1:size(i0,3)
        %A(:,:,i) = A(:,:,i)/magFactor(i)^2;
        A(:,:,i) = A(:,:,i)./tmpFac(:,:,i);
    end
    if bt
        specProb = specData.tweaker(:,:,k);
        broadTmp = specProb.*A.*i0.*exp(-atten*x2).*(eDen.^x3);
    else
        broadTmp = specProb(k)*A.*i0.*exp(-atten*x2).*(eDen.^x3);
    end
    broad = broad+broadTmp;
end

if (gamma)
    broad = edge_factor(scatFac.eFac.*gamma(1),eDen,broad,abs(gamma(2)),cg);
end

gBroad = zeros(size(scat1));
for i = 1:size(i0,3)
    tmpC = sqrt(magFactor(i))*35;
    gBroad(:,:,i) = exp(-(us.^2+vs.^2)/(tmpC^2));
end
if collim
    fg = fft2(gBroad.*triang);
else
    fg = fft2(gBroad);
end
scat1 = scat1+fft2(padarray(broad,[nPad,nPad,0])).*fg;

for j = 1:size(projA,3)
        projScat(:,:,j) = real(unpad(ifftshift(ifft2(scat1(:,:,j))),nPad));
end

projScat(projScat<0) = 0;
% n = projScat;
% b = projScat2;
out = projScat;
end

function out = unpad(in,nPad)
sz1 = size(in,1); sz2 = size(in,2);
out = in(1+nPad:sz1-nPad,1+nPad:sz2-nPad,:);
end

function out = edge_factor(factor,eDen,broadIn,sig,cg)
eDen = imgaussfilt(eDen,10);
[gradX,gradY] = gradient_op3d(eDen);

shiftFacX = 0.1552*eDen.*gradX.*factor./cg.ds;
shiftFacY = 0.1552*eDen.*gradY.*factor./cg.dt;
shiftFacX(abs(shiftFacX)>sig) = sig;
shiftFacY(abs(shiftFacX)>sig) = sig;

out = broadIn.*exp(-(shiftFacX.^2)./(35.^2)-(shiftFacY.^2)./(35.^2));
end
