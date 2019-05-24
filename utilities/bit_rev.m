function out = bit_rev(in,sz)
facs = factor(sz);
unit = zeros(size(facs));
unit(1) = 1;
if length(facs)>1
    for i = 2:length(facs)
        unit(i) = prod(facs(1:i-1));
    end
end
remain = in;
vec = zeros(size(facs));
while remain>0
    ind = 1;
    while unit(ind)*facs(ind)-1<remain
        ind = ind+1;
    end
    vec(ind) = vec(ind)+1;
    remain = remain-unit(ind);
end
out = vec(end);
ind = 1;
while ind<length(facs)
    out = out+prod(facs(end-ind+1:end))*vec(end-ind);
    ind = ind+1;
end
end