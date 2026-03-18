import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import label
from .sdf_superquadric import eul2rotm, rotm2eul


def _parse_input_args(grid, **kwargs):
    """Parse and set default hyperparameters for MPS algorithm."""
    defaults = {
        'verbose': True,
        'paddingSize': int(np.ceil(12 * grid['truncation'] / grid['interval'])),
        'minArea': int(np.ceil(grid['size'][0] / 20)),
        'maxDivision': 50,
        'scaleInitRatio': 0.1,
        'nanRange': 0.5 * grid['interval'],
        'w': 0.99,
        'tolerance': 1e-6,
        'relative_tolerance': 1e-4,
        'switch_tolerance': 1e-1,
        'maxSwitch': 2,
        'iter_min': 2,
        'maxOptiIter': 3,
        'maxIter': 15,
        'activeMultiplier': 3,
    }
    defaults.update(kwargs)
    return defaults


def _idx2coordinate(idx, grid):
    """
    Convert grid indices to 3D coordinates.

    Parameters
    ----------
    idx : ndarray (3, N)
        Grid indices (1-based, as in MATLAB).
    grid : dict
        Grid metadata.

    Returns
    -------
    coordinate : ndarray (3, N)
    """
    idx_floor = np.floor(idx).astype(int)
    idx_floor[idx_floor == 0] = 1

    # Convert to 0-based for array indexing
    x = grid['x'][idx_floor[0] - 1] + (idx[0] - idx_floor[0]) * grid['interval']
    y = grid['y'][idx_floor[1] - 1] + (idx[1] - idx_floor[1]) * grid['interval']
    z = grid['z'][idx_floor[2] - 1] + (idx[2] - idx_floor[2]) * grid['interval']
    return np.array([x, y, z])


def _idx3d_flatten(idx3d, grid):
    """
    Flatten 3D indices (1-based) to 1D indices (0-based for Python).

    Parameters
    ----------
    idx3d : ndarray (3, N)
        3D indices [row, col, slice] (1-based).
    grid : dict

    Returns
    -------
    idx : ndarray (N,)
        Flattened 1D indices (0-based).
    """
    sz = grid['size']
    # MATLAB: idx = x + size(1)*(y-1) + size(1)*size(2)*(z-1) => 1-based
    # Python: 0-based
    return ((idx3d[0] - 1) + sz[0] * (idx3d[1] - 1) +
            sz[0] * sz[1] * (idx3d[2] - 1)).astype(int)


def _regionprops(binary_mask):
    """
    Compute region properties similar to MATLAB regionprops.

    Parameters
    ----------
    binary_mask : ndarray (3D bool)

    Returns
    -------
    regions : list of dict
        Each dict has 'PixelIdxList', 'Area', 'Centroid', 'BoundingBox'.
    """
    # Use 26-connectivity (full) to match MATLAB's regionprops default for 3D
    struct_26 = np.ones((3, 3, 3), dtype=int)
    labeled, num_features = label(binary_mask, structure=struct_26)
    regions = []

    for i in range(1, num_features + 1):
        # Get 3D coordinates of the region
        coords = np.array(np.where(labeled == i))  # (3, N) in [row, col, slice]
        area = coords.shape[1]

        # F-order flat indices (matching SDF flat array layout from MATLAB)
        pixel_indices = np.ravel_multi_index(coords, binary_mask.shape, order='F')
        # Centroid (1-based, MATLAB convention: [col, row, slice])
        centroid = np.mean(coords, axis=1) + 1  # [row, col, slice] 0-based -> 1-based
        centroid_matlab = np.array([centroid[1], centroid[0], centroid[2]])  # [col, row, slice]

        # Bounding box (MATLAB convention: [col_start, row_start, slice_start, width, height, depth])
        mins = np.min(coords, axis=1)  # [row, col, slice] 0-based
        maxs = np.max(coords, axis=1)
        bbox = np.array([
            mins[1] + 1,  # col start (1-based)
            mins[0] + 1,  # row start (1-based)
            mins[2] + 1,  # slice start (1-based)
            maxs[1] - mins[1] + 1,  # width
            maxs[0] - mins[0] + 1,  # height
            maxs[2] - mins[2] + 1,  # depth
        ], dtype=float)

        regions.append({
            'PixelIdxList': pixel_indices,
            'Area': area,
            'Centroid': centroid_matlab,
            'BoundingBox': bbox,
        })

    return regions


def _rotz(deg):
    """Rotation matrix about Z axis, angle in degrees."""
    x = deg / 180.0 * np.pi
    return np.array([
        [np.cos(x), -np.sin(x), 0],
        [np.sin(x),  np.cos(x), 0],
        [0, 0, 1]
    ])


def _safe_power(base, exp):
    """Compute base**exp handling numerical edge cases like MATLAB."""
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        result = np.power(base, exp)
    return np.nan_to_num(result, nan=0.0, posinf=1e30, neginf=-1e30)


def _compute_sq_sdf(para, points):
    """Core superquadric SDF formula shared by internal and residual functions."""
    R = eul2rotm(para[5:8])
    t = para[8:11]
    X = R.T @ points - (R.T @ t[:, np.newaxis])

    r0 = np.linalg.norm(X, axis=0)

    eps1, eps2 = max(para[0], 1e-10), max(para[1], 1e-10)
    ax, ay, az = max(para[2], 1e-10), max(para[3], 1e-10), max(para[4], 1e-10)

    term_xy = _safe_power((np.abs(X[0]) / ax) ** 2, 1.0 / eps2) + \
              _safe_power((np.abs(X[1]) / ay) ** 2, 1.0 / eps2)
    term_z = _safe_power((np.abs(X[2]) / az) ** 2, 1.0 / eps1)
    scale = _safe_power(
        _safe_power(term_xy, eps2 / eps1) + term_z,
        -eps1 / 2.0
    )

    sdf = r0 * (1.0 - scale)
    sdf = np.nan_to_num(sdf, nan=0.0, posinf=1e30, neginf=-1e30)
    return sdf


def _sdf_superquadric_internal(para, points, truncation):
    """Internal SDF computation for a superquadric (same as in sdf_superquadric.py)."""
    sdf = _compute_sq_sdf(para, points)
    if truncation != 0:
        sdf = np.clip(sdf, -truncation, truncation)
    return sdf


def _difference_sq_sdf(para, sdf, points, truncation, weight):
    """Compute weighted SDF residual for least-squares fitting."""
    sdf_para = _compute_sq_sdf(para, points)

    if truncation != 0:
        sdf_para = np.clip(sdf_para, -truncation, truncation)

    dist = (sdf_para - sdf) * np.sqrt(weight)
    return dist


def _inlier_weight(sdf_active, active_idx, sdf_current, sigma2, w, truncation):
    """Compute inlier weights for robust fitting."""
    in_idx = sdf_active < 0.0 * truncation
    sdf_current_active = sdf_current[active_idx]

    weight = np.ones_like(sdf_active)
    if sigma2 < 1e-30:
        return weight

    const = w / ((1 - w) * (2 * np.pi * sigma2) ** (-0.5) * 1 * truncation)
    dist_current = np.clip(sdf_current_active[in_idx], -truncation, truncation) - sdf_active[in_idx]

    p = np.exp(-1.0 / (2 * sigma2) * dist_current ** 2)
    p = p / (const + p)
    weight[in_idx] = p
    return weight


def _cost_switched(candidates, sdf, points, truncation, weight):
    """Compute cost for candidate configurations."""
    values = np.zeros(candidates.shape[0])
    for i in range(candidates.shape[0]):
        residual = _difference_sq_sdf(candidates[i], sdf, points, truncation, weight)
        values[i] = np.sum(residual ** 2)
    return values


def _fit_superquadric_tsdf(sdf, x_init, truncation, points, roi_idx,
                           bounding_points, para):
    """
    Fit a superquadric to a local SDF region.

    Returns
    -------
    x : ndarray (11,) - fitted parameters
    occ_idx : ndarray - occupied voxel indices
    valid : ndarray (6,) - validity flags
    num_idx : ndarray (3,) - count statistics
    """
    # Positional bounds
    t_lb = bounding_points[:, 0]
    t_ub = bounding_points[:, 7]

    lb = np.array([0.0, 0.0, truncation, truncation, truncation,
                   -2 * np.pi, -2 * np.pi, -2 * np.pi,
                   t_lb[0], t_lb[1], t_lb[2]])
    ub = np.array([2.0, 2.0, 1.0, 1.0, 1.0,
                   2 * np.pi, 2 * np.pi, 2 * np.pi,
                   t_ub[0], t_ub[1], t_ub[2]])

    x = x_init.copy()
    cost = 0.0
    switched = 0
    nan_idx = ~np.isnan(sdf)
    sigma2 = np.exp(truncation) ** 2
    valid = np.zeros(6)

    for iteration in range(1, para['maxIter'] + 1):
        Rot = eul2rotm(x[5:8])
        check_points = np.array([
            x[8:11] - Rot[:, 0] * x[2],
            x[8:11] + Rot[:, 0] * x[2],
            x[8:11] - Rot[:, 1] * x[3],
            x[8:11] + Rot[:, 1] * x[3],
            x[8:11] - Rot[:, 2] * x[4],
            x[8:11] + Rot[:, 2] * x[4]
        ])  # (6, 3)

        valid[:3] = np.min(check_points, axis=0) >= t_lb - truncation
        valid[3:6] = np.max(check_points, axis=0) <= t_ub + truncation

        if not np.all(valid):
            break

        # SDF of voxels to current superquadric
        sdf_current = _sdf_superquadric_internal(x, points, 0)

        active_idx = ((sdf_current < para['activeMultiplier'] * truncation) &
                      (sdf_current > -para['activeMultiplier'] * truncation) &
                      nan_idx)

        points_active = points[:, active_idx]
        sdf_active = sdf[active_idx]

        weight = _inlier_weight(sdf_active, active_idx, sdf_current, sigma2, para['w'], truncation)

        # Scale upper bound
        Rot = eul2rotm(x[5:8])
        bP = bounding_points - x[8:11, np.newaxis]
        bP_body = Rot.T @ bP
        scale_limit = np.mean(np.abs(bP_body), axis=1)
        ub[2:5] = scale_limit

        # Least-squares optimization
        def cost_func(p):
            return _difference_sq_sdf(p, sdf_active, points_active, truncation, weight)

        # Clip initial guess to bounds (MATLAB's lsqnonlin does this automatically)
        x_clipped = np.clip(x, lb, ub)
        result = least_squares(cost_func, x_clipped, bounds=(lb, ub),
                               method='trf', max_nfev=3)
        x_n = result.x
        cost_n = result.cost * 2  # least_squares returns 0.5*sum(r^2)

        # Update sigma
        n_active = len(sdf_active)
        sigma2_n = cost_n / max(np.sum(weight), 1e-10)

        # Average cost
        cost_n_avg = cost_n / max(n_active, 1)

        # Relative cost decrease
        relative_cost = abs(cost - cost_n_avg) / max(cost_n_avg, 1e-10)

        if (cost_n_avg < para['tolerance'] and iteration > 1) or \
           (relative_cost < para['relative_tolerance'] and
            switched >= para['maxSwitch'] and iteration > para['iter_min']):
            x = x_n
            break

        if relative_cost < para['switch_tolerance'] and iteration != 1 \
                and switched < para['maxSwitch']:
            # Switching algorithm to avoid local minima
            switch_success = False

            # Case 1: axis-mismatch similarity
            axis_0 = eul2rotm(x[5:8])
            axis_1 = np.roll(axis_0, 2, axis=1)
            axis_2 = np.roll(axis_0, 1, axis=1)
            eul_1 = rotm2eul(axis_1)
            eul_2 = rotm2eul(axis_2)

            x_axis = np.array([
                [x[1], x[0], x[3], x[4], x[2], *eul_1, *x[8:11]],
                [x[1], x[0], x[4], x[2], x[3], *eul_2, *x[8:11]]
            ])

            # Case 2: duality similarity
            scale_ratio = np.roll(x[2:5], 2) / x[2:5]
            scale_idx = np.where((scale_ratio > 0.8) & (scale_ratio < 1.2))[0]
            x_rot_list = []

            # scale_idx uses 0-based indexing; MATLAB uses 1-based
            if 0 in scale_idx:  # MATLAB: ismember(1, scale_idx)
                eul_rot = rotm2eul(axis_0 @ _rotz(45))
                if x[1] <= 1:
                    new_scale = ((1 - np.sqrt(2)) * x[1] + np.sqrt(2)) * min(x[2], x[3])
                else:
                    new_scale = ((np.sqrt(2) / 2 - 1) * x[1] + 2 - np.sqrt(2) / 2) * min(x[2], x[3])
                x_rot_list.append([x[0], 2 - x[1], new_scale, new_scale, x[4],
                                   *eul_rot, *x[8:11]])

            if 1 in scale_idx:  # MATLAB: ismember(2, scale_idx)
                eul_rot = rotm2eul(axis_1 @ _rotz(45))
                if x[0] <= 1:
                    new_scale = ((1 - np.sqrt(2)) * x[0] + np.sqrt(2)) * min(x[3], x[4])
                else:
                    new_scale = ((np.sqrt(2) / 2 - 1) * x[0] + 2 - np.sqrt(2) / 2) * min(x[3], x[4])
                x_rot_list.append([x[1], 2 - x[0], new_scale, new_scale, x[2],
                                   *eul_rot, *x[8:11]])

            if 2 in scale_idx:  # MATLAB: ismember(3, scale_idx)
                eul_rot = rotm2eul(axis_2 @ _rotz(45))
                if x[0] <= 1:
                    new_scale = ((1 - np.sqrt(2)) * x[0] + np.sqrt(2)) * min(x[4], x[2])
                else:
                    new_scale = ((np.sqrt(2) / 2 - 1) * x[0] + 2 - np.sqrt(2) / 2) * min(x[4], x[2])
                x_rot_list.append([x[1], 2 - x[0], new_scale, new_scale, x[3],
                                   *eul_rot, *x[8:11]])

            x_rot = np.array(x_rot_list) if x_rot_list else np.zeros((0, 11))
            x_candidate = np.vstack([x_axis, x_rot]) if x_rot.shape[0] > 0 else x_axis

            # Evaluate candidates
            cost_candidate = _cost_switched(
                x_candidate, sdf_active, points_active, truncation, weight)

            valid_mask = ~np.isnan(cost_candidate) & ~np.isinf(cost_candidate)
            cost_candidate = cost_candidate[valid_mask]
            x_candidate = x_candidate[valid_mask]

            if len(cost_candidate) > 0:
                sort_idx = np.argsort(cost_candidate)

                for i_cand in sort_idx:
                    # Scale upper bound
                    Rot_c = eul2rotm(x_candidate[i_cand, 5:8])
                    bP_c = bounding_points - x_candidate[i_cand, 8:11, np.newaxis]
                    bP_body_c = Rot_c.T @ bP_c
                    scale_limit_c = np.mean(np.abs(bP_body_c), axis=1)
                    ub[2:5] = scale_limit_c

                    x_cand_clipped = np.clip(x_candidate[i_cand], lb, ub)
                    result_sw = least_squares(
                        cost_func, x_cand_clipped, bounds=(lb, ub),
                        method='trf', max_nfev=3)
                    cost_switch = result_sw.cost * 2

                    if cost_switch / max(n_active, 1) < min(cost_n_avg, cost):
                        x = result_sw.x
                        cost = cost_switch / max(n_active, 1)
                        sigma2 = cost_switch / max(np.sum(weight), 1e-10)
                        switch_success = True
                        break

            if not switch_success:
                cost = cost_n_avg
                x = x_n
                sigma2 = sigma2_n
            switched += 1
        else:
            cost = cost_n_avg
            sigma2 = sigma2_n
            x = x_n

    # Final occupancy check
    sdf_occ = _sdf_superquadric_internal(x, points, 0)
    occ = sdf_occ < para['nanRange']
    occ_idx = roi_idx[occ]
    occ_in = sdf_occ <= 0

    num_idx = np.zeros(3)
    sdf_occ_in = sdf[occ_in]
    num_idx[0] = np.sum((sdf_occ_in <= 0) | np.isnan(sdf_occ_in))
    num_idx[1] = np.sum(sdf_occ_in > 0)
    num_idx[2] = np.sum(sdf_occ_in <= 0)

    # Final validity check
    Rot = eul2rotm(x[5:8])
    check_points = np.array([
        x[8:11] - Rot[:, 0] * x[2],
        x[8:11] + Rot[:, 0] * x[2],
        x[8:11] - Rot[:, 1] * x[3],
        x[8:11] + Rot[:, 1] * x[3],
        x[8:11] - Rot[:, 2] * x[4],
        x[8:11] + Rot[:, 2] * x[4]
    ])
    valid[:3] = np.min(check_points, axis=0) >= t_lb - truncation
    valid[3:6] = np.max(check_points, axis=0) <= t_ub + truncation

    return x, occ_idx, valid, num_idx


def MPS(sdf, grid, **kwargs):
    """
    Marching-Primitives Shape abstraction algorithm.

    Extracts superquadric-based shape abstractions from a signed distance field.

    Parameters
    ----------
    sdf : ndarray (N,)
        Flattened truncated signed distance field.
    grid : dict
        Grid metadata with keys: 'size', 'range', 'x', 'y', 'z',
        'points', 'interval', 'truncation'.
    **kwargs : optional
        Algorithm hyperparameters (verbose, paddingSize, minArea, etc.)

    Returns
    -------
    x : ndarray (K, 11)
        Superquadric parameters. Each row:
        [eps1, eps2, ax, ay, az, eul_z, eul_y, eul_x, tx, ty, tz]
    """
    para = _parse_input_args(grid, **kwargs)
    sdf = sdf.copy().astype(float)

    num_division = 1
    x_all = []

    dratio = 3.0 / 5.0
    conn_ratio = np.array([dratio ** i for i in range(9)])
    conn_ratio[0] = 1.0
    conn_pointer = 0

    num_region = 1
    sz = grid['size']

    while num_division < para['maxDivision']:
        # 1 - Connectivity Marching
        if conn_pointer != 0 and num_region != 0:
            conn_pointer = 0

        conn_threshold = conn_ratio[conn_pointer] * np.nanmin(sdf)
        if conn_threshold > -grid['truncation'] * 3e-1:
            break

        # Reshape SDF to 3D for connectivity analysis
        sdf3d_region = sdf.reshape(sz[0], sz[1], sz[2], order='F')

        # Region connectivity analysis
        roi_list = _regionprops(sdf3d_region <= conn_threshold)

        # Filter by minimum area
        roi_list = [r for r in roi_list if r['Area'] >= para['minArea']]
        num_region = len(roi_list)

        if para['verbose']:
            print(f"Number of regions: {num_region}")

        if num_region == 0:
            if conn_pointer < len(conn_ratio) - 1:
                conn_pointer += 1
                continue
            else:
                break

        # 2 - Probabilistic Primitive Marching
        x_temp = np.zeros((num_region, 11))
        del_idx = np.zeros(num_region, dtype=bool)
        occ_idx_in = [None] * num_region
        num_idx_arr = np.zeros((num_region, 3))

        for i in range(num_region):
            roi = roi_list[i]
            occ_idx = np.array([], dtype=int)

            # Padding bounding box
            idx = np.ceil(roi['BoundingBox']).astype(int)
            padding = para['paddingSize']

            idx_end = np.minimum(
                idx[:3] + idx[3:6] + padding,
                [sz[1], sz[0], sz[2]]
            )
            idx_start = np.maximum(idx[:3] - padding, 1)

            # Generate 3D index grid
            range_row = np.arange(idx_start[1], idx_end[1] + 1)
            range_col = np.arange(idx_start[0], idx_end[0] + 1)
            range_slc = np.arange(idx_start[2], idx_end[2] + 1)
            idx_x, idx_y, idx_z = np.meshgrid(range_row, range_col, range_slc, indexing='ij')
            indices = np.array([idx_x.ravel(), idx_y.ravel(), idx_z.ravel()])

            roi_idx = _idx3d_flatten(indices, grid)

            # Bounding corner points
            bp_indices = np.array([
                [idx_start[1], idx_start[1], idx_end[1], idx_end[1],
                 idx_start[1], idx_start[1], idx_end[1], idx_end[1]],
                [idx_start[0], idx_start[0], idx_start[0], idx_start[0],
                 idx_end[0], idx_end[0], idx_end[0], idx_end[0]],
                [idx_start[2], idx_end[2], idx_start[2], idx_end[2],
                 idx_start[2], idx_end[2], idx_start[2], idx_end[2]]
            ], dtype=float)
            bounding_points = _idx2coordinate(bp_indices, grid)

            # Centroid
            centroid = np.maximum(np.floor(roi['Centroid']).astype(int), 1)
            centroid_flatten = _idx3d_flatten(
                np.array([[centroid[1]], [centroid[0]], [centroid[2]]]), grid
            )[0]

            if centroid_flatten in roi['PixelIdxList']:
                centroid_coord = grid['points'][:, centroid_flatten]
            else:
                # Search for nearest inside point
                pixel_points = grid['points'][:, roi['PixelIdxList']].T
                centroid_pt = grid['points'][:, centroid_flatten]
                dists = np.linalg.norm(pixel_points - centroid_pt, axis=1)
                k = np.argmin(dists)
                centroid_coord = grid['points'][:, roi['PixelIdxList'][k]]

            bbox_for_scale = np.array([
                idx_start[0], idx_start[1], idx_start[2],
                idx_end[0], idx_end[1], idx_end[2]
            ], dtype=float)

            valid = np.zeros(6)
            while not np.all(valid):
                # Initialize superquadric scale
                scale_init = para['scaleInitRatio'] * \
                    (bbox_for_scale[3:6] - bbox_for_scale[0:3]) * grid['interval']

                x_init = np.array([
                    1.0, 1.0,
                    scale_init[1], scale_init[0], scale_init[2],
                    0.0, 0.0, 0.0,
                    centroid_coord[0], centroid_coord[1], centroid_coord[2]
                ])

                # Fit superquadric
                x_temp[i], occ_idx, valid, num_idx_arr[i] = _fit_superquadric_tsdf(
                    sdf[roi_idx], x_init, grid['truncation'],
                    grid['points'][:, roi_idx], roi_idx,
                    bounding_points, para
                )

                if not np.all(valid):
                    extense = ~valid.astype(bool)
                    # Swap indices like MATLAB
                    extense_swapped = extense.copy()
                    extense_swapped[0], extense_swapped[1] = extense[1], extense[0]
                    extense_swapped[3], extense_swapped[4] = extense[4], extense[3]

                    # Check if at boundary
                    at_boundary = False
                    if any([
                        any(idx_start[extense_swapped[:3]] == 1),
                        idx_end[0] == sz[1] and extense_swapped[3],
                        idx_end[1] == sz[0] and extense_swapped[4],
                        idx_end[2] == sz[2] and extense_swapped[5]
                    ]):
                        at_boundary = True

                    if at_boundary:
                        break

                    # Extend bounding box
                    idx_extend = (~valid.astype(bool)).astype(int) * padding
                    idx_end = np.minimum(
                        idx_end + np.array([idx_extend[4], idx_extend[3], idx_extend[5]]),
                        [sz[1], sz[0], sz[2]]
                    )
                    idx_start = np.maximum(
                        idx_start - np.array([idx_extend[1], idx_extend[0], idx_extend[2]]),
                        1
                    )

                    # Recompute indices
                    range_row = np.arange(idx_start[1], idx_end[1] + 1)
                    range_col = np.arange(idx_start[0], idx_end[0] + 1)
                    range_slc = np.arange(idx_start[2], idx_end[2] + 1)
                    idx_x, idx_y, idx_z = np.meshgrid(
                        range_row, range_col, range_slc, indexing='ij')
                    indices = np.array([idx_x.ravel(), idx_y.ravel(), idx_z.ravel()])
                    roi_idx = _idx3d_flatten(indices, grid)

                    bp_indices = np.array([
                        [idx_start[1], idx_start[1], idx_end[1], idx_end[1],
                         idx_start[1], idx_start[1], idx_end[1], idx_end[1]],
                        [idx_start[0], idx_start[0], idx_start[0], idx_start[0],
                         idx_end[0], idx_end[0], idx_end[0], idx_end[0]],
                        [idx_start[2], idx_end[2], idx_start[2], idx_end[2],
                         idx_start[2], idx_end[2], idx_start[2], idx_end[2]]
                    ], dtype=float)
                    bounding_points = _idx2coordinate(bp_indices, grid)

                    bbox_for_scale = np.array([
                        idx_start[0], idx_start[1], idx_start[2],
                        idx_end[0], idx_end[1], idx_end[2]
                    ], dtype=float)

            # Store inside occupied indices
            occ_in_mask = sdf[occ_idx] <= 0 if len(occ_idx) > 0 else np.array([], dtype=bool)
            occ_idx_in[i] = occ_idx[occ_in_mask] if len(occ_idx) > 0 else np.array([], dtype=int)

        # Evaluate fitting quality
        for i in range(num_region):
            roi = roi_list[i]
            total = num_idx_arr[i, 0] + num_idx_arr[i, 1]
            out_pct = num_idx_arr[i, 1] / total if total > 0 else 1.0

            if out_pct > 0.3 or num_idx_arr[i, 0] < para['minArea'] or num_idx_arr[i, 2] <= 1:
                del_idx[i] = True
                sdf[roi['PixelIdxList']] = np.nan
                if para['verbose']:
                    print(f"region {i + 1}/{num_region}"
                          f" outPercentage: {out_pct:.4f}"
                          f" inNumber: {num_idx_arr[i, 2]:.0f}"
                          f" ...REJECTED")
            else:
                if occ_idx_in[i] is not None and len(occ_idx_in[i]) > 0:
                    sdf[occ_idx_in[i]] = np.nan
                if para['verbose']:
                    print(f"region {i + 1}/{num_region}"
                          f" outPercentage: {out_pct:.4f}"
                          f" inNumber: {num_idx_arr[i, 2]:.0f}"
                          f" ...ACCEPTED")

        # Collect accepted results
        accepted = x_temp[~del_idx]
        if len(accepted) > 0:
            x_all.append(accepted)
        num_division += 1

    if x_all:
        return np.vstack(x_all)
    else:
        return np.zeros((0, 11))
