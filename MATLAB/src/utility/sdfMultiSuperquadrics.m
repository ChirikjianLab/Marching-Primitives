function [sdf_u] = sdfMultiSuperquadrics(x, points, truncation)

n = size(x, 1);
sdf_u = sdfSuperquadric(x(1, :), points, truncation);
if n > 1
    for i = 2 : n
        sdf_u = min(sdf_u, ...
            sdfSuperquadric(x(i, :), points, truncation));
        
    end
end

end