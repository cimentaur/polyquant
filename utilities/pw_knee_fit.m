function [out,fac] = pw_knee_fit(specData,attenuationDb,mden,tiMa)
comp = attenuationDb.comp(2:53,:);
comp = [comp,zeros(52,1)];
comp = [comp;zeros(1,14)];
comp(end,end) = 1;
dens = attenuationDb.density(2:53);
dens(end+1) = 4.506;
mA = attenuationDb.massAtten(2:53,:);
mA = [mA;tiMa'];
Z = [1,6,7,8,11,12,15,16,17,19,20,26,53,22];
A = [1.01,12.01,14.01,16.00,22.99,24.305,30.97,32.066,35.45,39.098,40.08,55.845,126.90,47.867];
eDen = comp*(Z'./A');
waterE = 0.1119*1/1.01+0.8881*8/16;
[relDen,unIdx] = unique(dens.*eDen./waterE);
dens = dens(unIdx);
monoAtten = mA(unIdx,9).*dens;
if mden, relDen = dens; end
if mden == 2, relDen = monoAtten; end
warning('off','MATLAB:rankDeficientMatrix');
obj = @(z) ls_fit(specData,mA(unIdx,:),dens,relDen,z);
if mden == 2
    fac = fminsearch(obj,[0.15;1]);
    [~,out] = ls_fit(specData,mA(unIdx,:),dens,relDen,fac);
    %fac = [0.2;1];
else
    fac = fminsearch(obj,[1;3]);
    [~,out] = ls_fit(specData,mA(unIdx,:),dens,relDen,fac);
end
end

function [out,knee] = ls_fit(specData,mA,dens,relDen,z)
fac = [z;5];
res = zeros(size(dens));
for k = 1:length(specData.energy)
    fit = ls_knee(relDen,mA(:,k).*dens,fac);
    intPoint = interp1([0;fac],[0;fit],relDen);
    res = res+(mA(:,k).*dens-intPoint).^2;
    knee(:,k) = fac_to_knee([0;fac],[0;fit]);
end
res = sum(res);
out = res;
end

function out = ls_knee(in,y,kneeArray)
A = zeros(length(in),length(kneeArray)+1);
yOrd = zeros(size(y));
ind = [0;kneeArray];
start = 1;
for k = 1:length(kneeArray)
    inSub = in(in<ind(k+1) & in>=ind(k));
    ySub = y(in<ind(k+1) & in>=ind(k));
    aSub1 = (-inSub+ind(k))./(ind(k+1)-ind(k))+1;
    aSub2 = (inSub-ind(k))./(ind(k+1)-ind(k));
    A(start:(start+length(inSub)-1),[k,k+1]) = [aSub1,aSub2];
    yOrd(start:(start+length(inSub)-1)) = ySub;
    start = start+length(inSub);
end
A(:,1) = [];
out = A\yOrd;
end

function out = fac_to_knee(fac,fit)
knee = zeros((length(fac)-1)*2,1);
for k = 1:length(fit)-1
    knee(2*(k-1)+1) = (fit(k+1)-fit(k))/(fac(k+1)-fac(k));
    if k>1
        knee(2*k) = fit(k)-fac(k)*(knee(2*(k-1)+1));
    else
        knee(2*k) = 0;
    end
end
knee(2) = [];
out = knee;
end