function [out,flag] = ell_centroid(theta,rx,ry,sx,sy,y,flagIt)
m = y./150;
c = y-50*m;
rx = rx.^2; ry = ry.^2;
A = (cosd(theta)-m*sind(theta)).^2./rx+(sind(theta)+m*cosd(theta)).^2./ry;
B = -2*(c*sind(theta)+sx).*(cosd(theta)-m*sind(theta))./rx+...
    2*(c*cosd(theta)-sy).*(sind(theta)+m*cosd(theta))./ry;
C = (c*sind(theta)+sx).^2./rx+(c*cosd(theta)-sy).^2./ry-1;


x1 = (-B-sqrt(B.^2-4*A.*C))./(2*A);
x2 = (-B+sqrt(B.^2-4*A.*C))./(2*A);
flag = B.^2<4*A.*C;
if flagIt
    
    okInd = find(~flag);
    x1(1:min(okInd)) = x1(min(okInd));
    x1(max(okInd):end) = x1(max(okInd));
    x2(1:min(okInd)) = x2(min(okInd));
    x2(max(okInd):end) = x2(max(okInd));
end
out = (x1+x2)./2;
end