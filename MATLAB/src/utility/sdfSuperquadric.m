function [sdf] = sdfSuperquadric(para, points, truncation)

R = eul2rotm(para(6 : 8));
t = para(9 : 11);
X = R' * points - R' * t';

r0 = vecnorm(X);
scale = ((((X(1, :) / para(3)) .^ 2) .^ (1 / para(2)) + ...
        ((X(2, :) / para(4)) .^ 2) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
        ((X(3, :) / para(5)) .^ 2) .^ (1 / para(1))) .^ (-para(1) / 2);
        
sdf = r0 .* (1 - scale);

if truncation ~= 0
    sdf = min(max(sdf, - truncation), truncation);
end
end