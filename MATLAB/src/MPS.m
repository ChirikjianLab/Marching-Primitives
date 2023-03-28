function [x] = MPS(sdf, grid, varargin)
% Weixiao Liu 2023 Johns Hopkins University/National University of Singapore
% {Parameters--------------------------------------------------------------
% sdf is a vector containing the flattened truncated signed distance field
% grid is a object containing the gridding info of the TSDF
% para containes the algorithm hyper parameters

% Parsing varagin----------------------------------------------------------
[para] = parseInputArgs(grid, varargin{:});

% Initialization-----------------------------------------------------------
numDivision = 1; % initialize the depth of subdivision
x = []; % declare the array storing the superquadric configuration
dratio = 4/5; %4/5
connRatio = [1, dratio, (dratio)^2, (dratio)^3, (dratio)^4, ...
             (dratio)^5, (dratio)^6, (dratio)^7, (dratio)^8];
connPointer = 1;
num_region = 1;
occ = {};

while numDivision < para.maxDivision 
% terminate when the maximum layer of division reached 
% 1-Connectivity Marching--------------------------------------------------
    if connPointer ~= 1 && num_region ~= 0
        connPointer = 1;
    end
    % set the connection region 
    connTreshold = connRatio(connPointer) * min(sdf);
    if connTreshold > -grid.truncation * 1e-2
        break
    end
    % preliminary segementation via connectivity of regions
    % rearrange sdf to 3D array for region connection checking
    sdf3d_region = reshape(sdf, grid.size(1), grid.size(2), grid.size(3));
    % connection checking and preliminary subdivision
    roi = regionprops(sdf3d_region <= connTreshold, ...
        "PixelIdxList", "Area", "Centroid", "BoundingBox");
    
    % calculate the size of interest regions 
    roi = roi([roi.Area] >= para.minArea);    
    num_region = size(roi, 1);
    
    if para.verbose == 1
        disp(['Number of regions: ', num2str(num_region)])
    end
    
    if num_region == 0
        if connPointer < size(connRatio, 2)
            connPointer = connPointer + 1;
            continue
        else
            break
        end
    end           

% 2-Probabilistic Primitive Marching---------------------------------------  
    % initialize superquadric storage
    x_temp = zeros(num_region, 11);
    del_idx = zeros(num_region, 1);
    
    occ_temp = cell(num_region, 1);
    
    for i = 1 : num_region
        % padding the bounding box to extend the region of interest
        % also update the bounding box info accordingly
        idx = ceil(roi(i).BoundingBox);
        idx(4 : 6) = min(idx(1 : 3) + idx(4 : 6) + para.paddingSize, ...
            [grid.size(2), grid.size(1), grid.size(3)]);
        idx(1 : 3) = max(idx(1 : 3) - para.paddingSize, 1);
        [idx_x, idx_y, idx_z] = ndgrid( ...
            idx(2) : idx(5), idx(1) : idx(4), idx(3) : idx(6));
        roi(i).BoundingBox = idx;
        % find the voxels activated, that is, within the bounding box
        indices = reshape(cat(4, idx_x, idx_y, idx_z), [], 3)';
        % indices that is activated in each subdivision
        roi(i).idx = idx3d_flatten(indices, grid);
        roi(i).BoundingPoints = idx2Coordinate(...
            [idx(2), idx(2), idx(5), idx(5),...
             idx(2), idx(2), idx(5), idx(5);...
             idx(1), idx(1), idx(1), idx(1), ...
             idx(4), idx(4), idx(4), idx(4); ...
             idx(3), idx(6), idx(3), idx(6),...
             idx(3), idx(6), idx(3), idx(6)],...
             grid);
        
        % centralize the initial position inside the volumetric region
        roi(i).Centroid = max(floor(roi(i).Centroid), 1);                            
        centroid_flatten = idx3d_flatten(...
            [roi(i).Centroid(2); roi(i).Centroid(1); ...
            roi(i).Centroid(3)], grid);
        
        if ismember(centroid_flatten, roi(i).PixelIdxList)
            roi(i).Centroid = grid.points(:, centroid_flatten);
        else
            % search for nearest inside point within the roi
            k = dsearchn(grid.points(:, roi(i).PixelIdxList)', ...
                         grid.points(:, centroid_flatten)');
            roi(i).Centroid = grid.points(:, roi(i).PixelIdxList(k));
        end                           
              
        % initialize superquadric scale
        scale_init = para.scaleInitRatio ...
            * (roi(i).BoundingBox(4 : 6) - ...
            roi(i).BoundingBox(1 : 3)) * grid.interval;
        
        x_init = [1, 1, scale_init([2,1,3]) ...
            0, 0, 0, roi(i).Centroid'];
        
        valid = zeros(1 : 6);
        while ~all(valid)
            % for each subdivision find the optimal superquadric representation
            [x_temp(i, :), occ_idx, valid, num_idx] = fitSuperquadricTSDF( ...
                sdf(roi(i).idx), ...
                x_init, ...
                grid.truncation, ...
                grid.points(:, roi(i).idx), ...
                roi(i).idx, ...
                roi(i).BoundingPoints, ...
                para);

            if ~all(valid)
                
                extense = ~valid;
                extense([1, 2]) = extense([2, 1]);
                extense([4, 5]) = extense([5, 4]);
                if any([any(idx(extense(1 : 3)) == 1), ...
                       idx(4) == grid.size(2) && extense(4), ...
                       idx(5) == grid.size(1) && extense(5), ...
                       idx(6) == grid.size(3) && extense(6)])
                   break
                end               
                idx_extend = ~valid * para.paddingSize;
                idx(4 : 6) = min(idx(4 : 6) + idx_extend([5, 4, 6]), ...
                    [grid.size(2), grid.size(1), grid.size(3)]);
                idx(1 : 3) = max(idx(1 : 3) - idx_extend([2, 1, 3]), 1);
                [idx_x, idx_y, idx_z] = ndgrid( ...
                    idx(2) : idx(5), idx(1) : idx(4), idx(3) : idx(6));
                roi(i).BoundingBox = idx;
                % find the voxels activated, that is, within the bounding box
                indices = reshape(cat(4, idx_x, idx_y, idx_z), [], 3)';
                % indices that is activated in each subdivision
                roi(i).idx = idx3d_flatten(indices, grid);
                roi(i).BoundingPoints = idx2Coordinate(...
                    [idx(2), idx(2), idx(5), idx(5),...
                    idx(2), idx(2), idx(5), idx(5);...
                    idx(1), idx(1), idx(1), idx(1), ...
                    idx(4), idx(4), idx(4), idx(4); ...
                    idx(3), idx(6), idx(3), idx(6),...
                    idx(3), idx(6), idx(3), idx(6)],...
                    grid);                
            end
        end
        
        % evaluate fitting quality
        occ_idx_in = occ_idx(sdf(occ_idx) <= 0);       
        occ_temp{i} = occ_idx(sdf(occ_idx) < grid.truncation);
        
        if (num_idx(2) / (num_idx(1) + num_idx(2))) > 0.25 ... 
                || num_idx(3) < para.minArea %3
            del_idx(i) = 1;
            sdf(roi(i).PixelIdxList) = nan;
            if para.verbose == 1
                disp(['region ', num2str(i), '/', num2str(num_region), ...
                    ' outPrecentage: ',...
                    num2str(num_idx(2) / (num_idx(1) + num_idx(2))), ...
                    ' inNumber: ', ...
                    num2str(num_idx(3)), ...
                    ' ...REJECTED'])
            end
        else
            sdf(occ_idx_in) = nan;
            if para.verbose == 1
                disp(['region ', num2str(i), '/', num2str(num_region), ...
                    ' outPercentage: ',...
                    num2str(num_idx(2) / (num_idx(1) + num_idx(2))), ...
                    ' inNumber: ', ...
                    num2str(num_idx(3)), ...
                    '...ACCEPTED'])
            end
        end
    end

    x_temp = x_temp(del_idx == 0, :);
    x(end + 1 : end + size(x_temp, 1), :) = x_temp;
    numDivision = numDivision + 1;
    
    occ_temp = occ_temp(del_idx == 0, :);
    occ(end + 1 : end + size(occ_temp, 1), :) = occ_temp;
    
end
end

% Parse varargin-----------------------------------------------------------
function [para] = parseInputArgs(grid, varargin)

    % set input parser
    defaults = struct(...
        'verbose', true, ...
        'paddingSize', ...
        ceil(15 * grid.truncation / grid.interval),... 
        'minArea', 5, ...
        'maxDivision', 50, ...
        'scaleInitRatio', 0.1, ...
        'nanRange', 1.5 * grid.interval, ...
        'w', 0.99, ...
        'tolerance', 1e-6, ...
        'relative_tolerance', 1e-1, ...
        'switch_tolerance', 1e-1, ...
        'maxSwitch', uint8(2), ...
        'iter_min', uint8(3), ...
        'maxOptiIter', 2, ...
        'maxIter', 10, ...
        'activeMultiplier', 3.0); 

    parser = inputParser;
    parser.CaseSensitive = false;
    
    parser.addParameter('verbose', defaults.verbose, @islogical);
    parser.addParameter('paddingSize', defaults.paddingSize, ...
        @(x)validateattributes(x, {'single', 'double'}, {'integer', '>=', 0}));
    parser.addParameter('minArea', defaults.minArea, ...
        @(x)validateattributes(x, {'single', 'double'}, {'integer', '>=', 0}));
    parser.addParameter('maxDivision', defaults.maxDivision, ...
        @(x)validateattributes(x, {'single', 'double'}, {'integer', '>', 0}));
    parser.addParameter('scaleInitRatio', defaults.scaleInitRatio, ...
        @(x)validateattributes(x, {'single', 'double'}, {'>', 0}));
    parser.addParameter('nanRange', defaults.nanRange, ...
        @(x)validateattributes(x, {'single', 'double'}, {'>=', 0}));
    parser.addParameter('w', defaults.w, ...
        @(x)validateattributes(x, {'single', 'double'}, {'>=', 0, '<', 1}));
    parser.addParameter('tolerance', defaults.tolerance, ...
        @(x)validateattributes(x, {'single', 'double'}, {'>', 0}));
    parser.addParameter('relative_tolerance', defaults.relative_tolerance, ...
        @(x)validateattributes(x, {'single', 'double'}, {'>', 0}));
    parser.addParameter('switch_tolerance', defaults.switch_tolerance, ...
        @(x)validateattributes(x, {'single', 'double'}, {'>', 0}));
    parser.addParameter('maxSwitch', defaults.maxSwitch, ...
        @(x)validateattributes(x, {'single', 'double'}, {'integer', '>=', 0}));
    parser.addParameter('iter_min', defaults.iter_min, ...
        @(x)validateattributes(x, {'single', 'double'}, {'integer', '>=', 0}));
    parser.addParameter('maxOptiIter', defaults.maxOptiIter, ...
        @(x)validateattributes(x, {'single', 'double'}, {'integer', '>', 0}));
    parser.addParameter('maxIter', defaults.maxIter, ...
        @(x)validateattributes(x, {'single', 'double'}, {'integer', '>', 0}));
    parser.addParameter('activeMultiplier', defaults.activeMultiplier, ...
        @(x)validateattributes(x, {'single', 'double'}, {'>', 0}));    

    parser.parse(varargin{:});   
    para = parser.Results;   
end


% From grid index to 3D coordinates----------------------------------------
function [coordinate] = idx2Coordinate(idx, grid)
idx_floor = floor(idx);
idx_floor(idx_floor == 0) = 1;

x = grid.x(idx_floor(1, :)) + (idx(1, :) - idx_floor(1, :)) * grid.interval;
y = grid.y(idx_floor(2, :)) + (idx(2, :) - idx_floor(2, :)) * grid.interval;
z = grid.z(idx_floor(3, :)) + (idx(3, :) - idx_floor(3, :)) * grid.interval;
coordinate = [x; y; z];
end

% Flatten 3D index into 1D-------------------------------------------------
function [idx] = idx3d_flatten(idx3d, grid)
idx = reshape(idx3d(1, :) + grid.size(1) * (idx3d(2, :) - 1) + ...
    grid.size(1) * grid.size(2) * (idx3d(3, :) - 1), [], 1);
end

% Recover superquadric from SDF--------------------------------------------
function [x, occ_idx, valid, num_idx] = fitSuperquadricTSDF(...
    sdf, x_init, truncation, points, roi_idx, boundingPoints, para)
% grid is only used for debug visualization
options = optimoptions(...
    'lsqnonlin', 'Algorithm', 'trust-region-reflective', ...
    'Display', 'off', 'MaxIterations', para.maxOptiIter);

% initialize validity vector
valid = zeros(1, 6);

% positional upper and lower bound
t_lb = boundingPoints(:, 1)';
t_ub = boundingPoints(:, 8)';

% setting upper and lower bound
% lb = [0.0, 0.0, 0.0001, 0.0001, 0.0001, ...
%     -2 * pi, -2 * pi, -2 * pi, t_lb];
lb = [0.0, 0.0, truncation, truncation, truncation, ...
    -2 * pi, -2 * pi, -2 * pi, t_lb];
ub = [2, 2, 1, 1, 1, ...
    2 * pi, 2 * pi, 2 * pi, t_ub];

% initialization
x = x_init;
cost = 0;
switched = uint8(0);
nan_idx = ~isnan(sdf);
sigma2 = truncation;

for iter = 1 : para.maxIter    
    Rot = eul2rotm(x(6 : 8));
    checkPoints = [x(9 : 11) - Rot(:, 1)' * x(3);
                   x(9 : 11) + Rot(:, 1)' * x(3);
                   x(9 : 11) - Rot(:, 2)' * x(4);
                   x(9 : 11) + Rot(:, 2)' * x(4);
                   x(9 : 11) - Rot(:, 3)' * x(5);
                   x(9 : 11) + Rot(:, 3)' * x(5)];
    valid(1 : 3) = min(checkPoints) >= t_lb - 1 * truncation;
    valid(4 : 6) = max(checkPoints) <= t_ub + 1 * truncation;
    
    % break the iteration when validity is violated
    if ~all(valid)
        break
    end
    
    % calculating the signed distance of voxels to the current superquadric
    sdf_current = sdfSuperquadric(x, points, 0);
    % find the voxels activated by the current superquadric configuration
    active_idx = and(and(sdf_current < para.activeMultiplier* truncation, ...
                     sdf_current > -para.activeMultiplier* truncation), ...
                     nan_idx);
                 
    points_active = points(:, active_idx);
    sdf_active = sdf(:, active_idx);
    
    [weight] = inlierWeight(...
            sdf_active, active_idx, sdf_current, sigma2, para.w, truncation);
    
    % scale upper bound
    Rot = eul2rotm(x(6 : 8));
    bP = boundingPoints - x(9 : 11)';
    bP_body = Rot' * bP;
    scale_limit = mean(abs(bP_body'));   %max mean
    ub(3 : 5) = scale_limit;
    
    cost_func = @(para) differenceSQSDF(...
        para, sdf_active, points_active, truncation, weight);
    [x_n, cost_n] = lsqnonlin(cost_func, x, lb, ub, options);
      
    % update sigma 
    sigma2_n = cost_n / sum(weight);
    
    % average cost
    cost_n = cost_n / length(sdf_active);
    
    %evaluate relative cost decrease
    relative_cost = abs(cost - cost_n) / cost_n;
    
    if (cost_n < para.tolerance && iter > 1) || ...
            (relative_cost < para.relative_tolerance ...
            && switched >= para.maxSwitch && iter > para.iter_min)
        x = x_n;
        break
    end
    
    if relative_cost < para.switch_tolerance && iter ~= 1 ...
            && switched < para.maxSwitch
        
        % activate switching algorithm to avoid local minimum
        switch_success = 0;
        
        % case1 - axis-mismatch similarity
        axis_0 = eul2rotm(x(6 : 8));
        axis_1 = circshift(axis_0, [0, 2]);
        axis_2 = circshift(axis_0, [0, 1]);
        eul_1 = rotm2eul(axis_1);
        eul_2 = rotm2eul(axis_2);
        x_axis = [x(2), x(1), x(4), x(5), x(3), eul_1, x(9 : 11); ...
            x(2), x(1), x(5), x(3), x(4), eul_2, x(9 : 11)];
        
        % case2 - duality similarity and combinations
        scale_ratio = circshift(x(3 : 5), 2) ./ x(3 : 5);
        scale_idx = find(and(scale_ratio > 0.6, scale_ratio < 1.4));
        x_rot = zeros(size(scale_idx, 2), 11);
        rot_idx = 1;
        
        if ismember(1, scale_idx)
            eul_rot = rotm2eul(axis_0 * rotz(45));
            if x(2) <= 1
                x_rot(rot_idx, :) = [x(1), 2 - x(2), ...
                    ((1 - sqrt(2)) * x(2) + sqrt(2)) * min(x(3), x(4)) ...
                    * ones(1, 2), x(5), eul_rot, x(9 : 11)];
            else
                x_rot(rot_idx, :) = [x(1), 2 - x(2), ...
                    ((sqrt(2)/2 - 1) * x(2) + 2 - sqrt(2)/2) * ...
                    min(x(3), x(4)) * ones(1, 2), x(5), eul_rot, x(9 : 11)];
            end
            rot_idx = rot_idx + 1;
        end
        
        if ismember(2, scale_idx)
            eul_rot = rotm2eul(axis_1 * rotz(45));
            if x(1) <= 1
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ...
                    ((1 - sqrt(2)) * x(1) + sqrt(2)) * min(x(4), x(5)) ...
                    * ones(1, 2), x(3), eul_rot, x(9 : 11)];
            else
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ...
                    ((sqrt(2)/2 - 1) * x(1) + 2 - sqrt(2)/2) * ...
                    min(x(4), x(5)) * ones(1, 2), x(3), eul_rot, x(9 : 11)];
            end
            rot_idx = rot_idx + 1;
        end
        
        if ismember(3, scale_idx)
            eul_rot = rotm2eul(axis_2 * rotz(45));
            if x(1) <= 1
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ...
                    ((1 - sqrt(2)) * x(1) + sqrt(2)) * min(x(5), x(3)) ...
                    * ones(1, 2), x(4), eul_rot, x(9 : 11)];
            else
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ...
                    ((sqrt(2)/2 - 1) * x(1) + 2 - sqrt(2)/2) * ...
                    min(x(5), x(3)) * ones(1, 2), x(4), eul_rot, x(9 : 11)];
            end
        end
        
        % generate candidate configuration list with cost
        x_candidate = [x_axis; x_rot];
        cost_candidate = cost_switched(...
            x_candidate, sdf_active, points_active, truncation, weight);
        
        idx_nan = find(...
            and(~isnan(cost_candidate), ~isinf(cost_candidate)));
        cost_candidate = cost_candidate(idx_nan);
        x_candidate = x_candidate(idx_nan, :);
        
        [~, idx] = sort(cost_candidate);
        for i_candidate = 1 : size(idx, 1)
            
            % scale upper  bound
            Rot = eul2rotm(x_candidate(idx(i_candidate), 6 : 8));
            bP = boundingPoints - x_candidate(idx(i_candidate), 9:11)';
            bP_body = Rot' * bP;
            scale_limit = mean(abs(bP_body'));
            ub(3 : 5) = scale_limit;        
            
            [x_switch, cost_switch] = lsqnonlin(...
                cost_func, x_candidate(idx(i_candidate), :), ...
                lb, ub, options);
            
            if cost_switch / length(sdf_active) < min(cost_n, cost)
                x = x_switch;
                cost = cost_switch / length(sdf_active);
                sigma2 = cost_switch / sum(weight);
                switch_success = 1;
                break
            end
        end
        
        if switch_success == 0           
            cost = cost_n;
            x = x_n;
            sigma2 = sigma2_n;
        end
        switched = switched + 1;
    else
        cost = cost_n;
        sigma2 = sigma2_n;
        x = x_n;
    end
end

sdf_occ = sdfSuperquadric(x, points, 0);
occ = sdf_occ < para.nanRange;
occ_idx = roi_idx(occ);
occ_in = sdf_occ <= 0;

num_idx = zeros(1, 3);
num_idx(1) = sum(or(sdf(occ_in) <= 0, isnan(sdf(occ_in))));
num_idx(2) = sum(sdf(occ_in) > 0);
num_idx(3) = sum(sdf(occ_in) <= 0);

% final check size validity
Rot = eul2rotm(x(6 : 8));
checkPoints = [x(9 : 11) - Rot(:, 1)' * x(3);
               x(9 : 11) + Rot(:, 1)' * x(3);
               x(9 : 11) - Rot(:, 2)' * x(4);
               x(9 : 11) + Rot(:, 2)' * x(4);
               x(9 : 11) - Rot(:, 3)' * x(5);
               x(9 : 11) + Rot(:, 3)' * x(5)];

valid(1 : 3) = min(checkPoints) >= t_lb - 1 * truncation;
valid(4 : 6) = max(checkPoints) <= t_ub + 1 * truncation;
%--------------------------------------------------------------------------
    function [R] = eul2rotm(eul)
        R = zeros(3,3,size(eul,1),'like',eul);
        ct = cos(eul);
        st = sin(eul);
        R(1,1,:) = ct(:,2).*ct(:,1);
        R(1,2,:) = st(:,3).*st(:,2).*ct(:,1) - ct(:,3).*st(:,1);
        R(1,3,:) = ct(:,3).*st(:,2).*ct(:,1) + st(:,3).*st(:,1);
        R(2,1,:) = ct(:,2).*st(:,1);
        R(2,2,:) = st(:,3).*st(:,2).*st(:,1) + ct(:,3).*ct(:,1);
        R(2,3,:) = ct(:,3).*st(:,2).*st(:,1) - st(:,3).*ct(:,1);
        R(3,1,:) = -st(:,2);
        R(3,2,:) = st(:,3).*ct(:,2);
        R(3,3,:) = ct(:,3).*ct(:,2);                   
    end
%--------------------------------------------------------------------------
    function [weight] = inlierWeight(...
            sdf_active, active_idx, sdf_current, sigma2, w, truncation)
        
        inIdx = sdf_active < 0.0 * truncation;
        sdf_current = sdf_current(active_idx);
        
        const = w / ((1 - w) * (2 * pi * sigma2) ^ (- 1 / 2) * ...
            1 * truncation);
        dist_current = min(max(sdf_current(inIdx), ...
            - truncation), truncation) - sdf_active(inIdx);
        
        weight = ones(size(sdf_active));
        p = exp(-1 / (2 * sigma2) * dist_current .^ 2);
        p = p ./ (const + p);
        weight(inIdx) = p;
    end
%--------------------------------------------------------------------------
    function [value] = cost_switched(para, sdf, points, truncation, weight)
        value = zeros(size(para, 1), 1);
        for i = 1 : size(para, 1)
            value(i, 1) = sum(differenceSQSDF(para(i, :), ...
                sdf, points, truncation, weight) .^ 2, 2);
        end
    end
%--------------------------------------------------------------------------
    function [dist] = differenceSQSDF(para, sdf, points, truncation, weight)
        
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11);
        X = R' * points - R' * t';
        
        r0 = vecnorm(X);
                
        scale = ((((X(1, :) / para(3)) .^ 2) .^ (1 / para(2)) + ...
            ((X(2, :) / para(4)) .^ 2) .^ ...
            (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X(3, :) / para(5)) .^ 2) .^ (1 / para(1))) .^ (-para(1) / 2);
        
        sdf_para = r0 .* (1 - scale);
        
        if truncation ~= 0
            sdf_para = min(max(sdf_para, -truncation), truncation);
        end
        
        dist = (sdf_para - sdf) .* weight .^ (1 / 2);     
    end
%--------------------------------------------------------------------------
    function [sdf] = sdfSuperquadric(para, points, truncation)
        
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11);
        X = R' * points - R' * t';
        
        r0 = vecnorm(X);
        scale = ((((X(1, :) / para(3)) .^ 2) .^ (1 / para(2)) + ...
            ((X(2, :) / para(4)) .^ 2) .^ ...
            (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X(3, :) / para(5)) .^ 2) .^ (1 / para(1))) .^ (-para(1) / 2);
        
        sdf = r0 .* (1 - scale);
        
        if truncation ~= 0
            sdf = min(max(sdf, -truncation), truncation);
        end
    end
%--------------------------------------------------------------------------

end