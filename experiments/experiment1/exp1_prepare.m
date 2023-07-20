save_file_path = 'experiments/experiment1/exp1_inputs.mat';

rng(1000, 'combRecursive');

marg_num = 10;

marg_weights = ones(marg_num, 1) / marg_num;

marg_vert_cell = cell(marg_num, 1);
marg_tri_cell = cell(marg_num, 1);
marg_dens_cell = cell(marg_num, 1);

marg_vert_cell{1} = ...
    [0, 20; ...
    -4, 12;
    4, 12];
marg_tri_cell{1} = ...
    [1, 2, 3];

marg_vert_cell{2} = ...
    [-6, 12; ...
    -2, 12; ...
    2, 12; ...
    6, 12; ...
    -2, 8; ...
    2, 8; ...
    -6, 4; ...
    -2, 4; ...
    2, 4; ...
    6, 4];
marg_tri_cell{2} = ...
    [1, 2, 5; ...
    2, 5, 6; ...
    2, 3, 6; ...
    3, 4, 6; ...
    7, 8, 5; ...
    8, 9, 5; ...
    9, 5, 6; ...
    10, 9, 6];

marg_vert_cell{3} = ...
    [-16, 16; ...
    -8, 16; ...
    -12, 12; ...
    -16, 8; ...
    -8, 8];
marg_tri_cell{3} = ...
    [1, 2, 3; ...
    1, 3, 4; ...
    2, 3, 5; ...
    3, 4, 5];

marg_vert_cell{4} = ...
    [-14, 5; ...
    -10, 5; ...
    -14, 1; ...
    -10, 1; ...
    -6, 1; ...
    -14, -3; ...
    -10, -3; ...
    -6, -3; ...
    -14, -7; ...
    -10, -7];
marg_tri_cell{4} = ...
    [1, 2, 4; ...
    1, 3, 4; ...
    2, 4, 5; ...
    3, 4, 7; ...
    3, 6, 7; ...
    4, 5, 8; ...
    4, 7, 8; ...
    6, 7, 10; ...
    6, 9, 10; ...
    7, 8, 10];

marg_vert_cell{5} = ...
    [-7, -9; ...
    -11, -11; ...
    -15, -13; ...
    -7, -13; ...
    -11, -15; ...
    -7, -17];
marg_tri_cell{5} = ...
    [1, 2, 4; ...
    2, 3, 5; ...
    2, 4, 5; ...
    4, 5, 6];

marg_vert_cell{6} = ...
    [-4, -11; ...
    4, -11; ...
    0, -13; ...
    -4, -15; ...
    4, -15];
marg_tri_cell{6} = ...
    [1, 3, 4; ...
    2, 3, 5];

marg_vert_cell{7} = ...
    [-2, -5; ...
    2, -5; ...
    -6, -9; ...
    -2, -9; ...
    2, -9; ...
    6, -9; ...
    0, -13];
marg_tri_cell{7} = ...
    [1, 3, 4; ...
    1, 2, 5; ...
    1, 4, 5; ...
    2, 5, 6; ...
    3, 4, 7; ...
    4, 5, 7; ...
    5, 6, 7];

marg_vert_cell{8} = ...
    [10, -5; ...
    6, -8; ...
    14, -8; ...
    10, -11; ...
    6, -14; ...
    14, -14; ...
    10, -17];
marg_tri_cell{8} = ...
    [1, 2, 4; ...
    1, 3, 4; ...
    4, 5, 7; ...
    4, 6, 7];

marg_vert_cell{9} = ...
    [10, 8; ...
    14, 10; ...
    10, 6; ...
    14, 6; ...
    10, 2; ...
    14, 2; ...
    10, -2; ...
    14, -2; ...
    10, -4; ...
    14, -6];
marg_tri_cell{9} = ...
    [1, 2, 4; ...
    1, 3, 4; ...
    3, 4, 6; ...
    3, 5, 6; ...
    5, 6, 8; ...
    5, 7, 8; ...
    7, 8, 9; ...
    8, 9, 10];

marg_vert_cell{10} = ...
    [12, 19; ...
    8, 17; ...
    16, 17; ...
    12, 15; ...
    8, 13; ...
    16, 13; ...
    12, 11];
marg_tri_cell{10} = ...
    [1, 2, 4; ...
    1, 3, 4; ...
    2, 4, 5; ...
    3, 4, 6; ...
    4, 5, 7; ...
    4, 6, 7];

for marg_id = 1:marg_num
    vert_num = size(marg_vert_cell{marg_id}, 1);
    marg_dens_cell{marg_id} = gamrnd(1, 2, vert_num, 1);
end

step_num = 7;
mesh_sizes = [0, 0; 2, 2; 1, 1; 2/3, 2/3; 2/4, 2/4; 2/5, 2/5; 2/6, 2/6];
plot_grid_x_num = 1000;
plot_grid_y_num = 1000;

plot_pt_x_cell = cell(marg_num, 1);
plot_pt_y_cell = cell(marg_num, 1);
plot_inside_cell = cell(marg_num, 1);
plot_pts_cell = cell(marg_num, 1);

testfuncs_cell = cell(step_num, 1);

testfunc_num_mat = zeros(step_num, marg_num);

for step_id = 1:step_num
    testfuncs_cell{step_id} = cell(marg_num, 1);
end

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_vertices = marg_vert_cell{marg_id};
    marg_triangles = marg_tri_cell{marg_id};
    marg_dens = marg_dens_cell{marg_id};
    marg_x_min = min(marg_vertices(:, 1));
    marg_x_max = max(marg_vertices(:, 1));
    marg_y_min = min(marg_vertices(:, 2));
    marg_y_max = max(marg_vertices(:, 2));

    marg = ProbMeas2D_CPWADens(marg_vertices, marg_triangles, marg_dens);
    marg_cell{marg_id} = marg;

    for step_id = 1:step_num
        marg_grid_pts_x = marg_x_min:mesh_sizes(step_id, 1):marg_x_max;
        marg_grid_pts_y = marg_y_min:mesh_sizes(step_id, 2):marg_y_max;
        [marg_grid_x, marg_grid_y] = meshgrid(marg_grid_pts_x, ...
            marg_grid_pts_y);
        marg_grid = [marg_vert_cell{marg_id}; ...
        [marg_grid_x(:), marg_grid_y(:)]];

        marg_grid = unique(marg_grid, 'stable', 'rows');
        inside_support = marg.checkIfInsideSupport(marg_grid);
        marg_grid = marg_grid(inside_support, :);

        region_num = size(marg_tri_cell{marg_id}, 1);
        inside_mat = marg.checkIfInsideTriangularRegion(marg_grid);
        marg_triangles_cell = cell(region_num, 1);

        for region_id = 1:region_num
            region_marg_indices = find(inside_mat(:, region_id));
            region_triangulation = delaunay(marg_grid( ...
                region_marg_indices, :));
            marg_triangles_cell{region_id} = ...
                region_marg_indices(region_triangulation);

            if size(marg_triangles_cell{region_id}, 2) == 1
                marg_triangles_cell{region_id} = ...
                    marg_triangles_cell{region_id}';
            end

            region_tri_num = size(marg_triangles_cell{region_id}, 1);
            keep_list = true(region_tri_num, 1);

            for tri_id = 1:region_tri_num
                cond_num = cond( ...
                    [marg_grid(marg_triangles_cell{region_id}(tri_id, ...
                    :), :)'; ones(1, 3)]);

                if cond_num > 1e12
                    keep_list(tri_id) = false;
                end
            end

            marg_triangles_cell{region_id} = ...
                marg_triangles_cell{region_id}(keep_list, :);
        end

        testfuncs_cell{step_id}{marg_id} = ...
            {marg_grid; marg_triangles_cell; false};

        testfunc_num_mat(step_id, marg_id) = size(marg_grid, 1);
    end

    plot_pt_x = linspace(marg_x_min, marg_x_max, plot_grid_x_num + 1);
    plot_pt_x = (plot_pt_x(1:end - 1) + plot_pt_x(2:end)) / 2;
    plot_pt_y = linspace(marg_y_min, marg_y_max, plot_grid_y_num + 1);
    plot_pt_y = (plot_pt_y(1:end - 1) + plot_pt_y(2:end)) / 2;

    plot_pt_x_cell{marg_id} = plot_pt_x;
    plot_pt_y_cell{marg_id} = plot_pt_y;

    [plot_grid_x, plot_grid_y] = meshgrid(plot_pt_x, plot_pt_y);

    plot_pts_cell{marg_id} = [plot_grid_x(:), plot_grid_y(:)];

    plot_inside_cell{marg_id} = marg.checkIfInsideSupport( ...
        plot_pts_cell{marg_id});
end

quality_poly_num = 3;

quality_poly_cell = cell(quality_poly_num, 1);

sq = [0, 0; 5, 0; 5, 1; 0, 1];
sq_x_min = min(sq(:, 1));
sq_x_max = max(sq(:, 1));
sq_y_min = min(sq(:, 2));
sq_y_max = max(sq(:, 2));
q_mesh_sizes = [1/2, 1/2; 1/3, 1/3; 1/4, 1/4; 1/5, 1/5; 1/6, 1/6; ...
    1/7, 1/7; 1/8, 1/8];

sq_grid_cell = cell(step_num, 1);
sq_tri_cell = cell(step_num, 1);

for step_id = 1:step_num
    [sq_grid_x, sq_grid_y] = meshgrid( ...
        sq_x_min:q_mesh_sizes(step_id, 1):sq_x_max, ...
        sq_y_min:q_mesh_sizes(step_id, 2):sq_y_max);
    sq_grid_cell{step_id} = [sq_grid_x(:), sq_grid_y(:)];
    sq_tri_cell{step_id} = delaunay(sq_grid_cell{step_id});
end

angles = [0; 2 * pi / 3;  4 * pi / 3];
shifts = [-2, 1; 2, 0.8; -0.2, -2.7];

quality_x_min = inf;
quality_x_max = -inf;
quality_y_min = inf;
quality_y_max = -inf;

for poly_id = 1:quality_poly_num
    angle = angles(poly_id);
    rot_mat = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    quality_poly_cell{poly_id} = sq * rot_mat + shifts(poly_id, :);

    quality_x_min = min(quality_x_min, ...
        min(quality_poly_cell{poly_id}(:, 1)));
    quality_x_max = max(quality_x_max, ...
        max(quality_poly_cell{poly_id}(:, 1)));
    quality_y_min = min(quality_y_min, ...
        min(quality_poly_cell{poly_id}(:, 2)));
    quality_y_max = max(quality_y_max, ...
        max(quality_poly_cell{poly_id}(:, 2)));
end


quality_testfuncs_cell = cell(step_num, 1);
quality_testfuncs_num_list = zeros(step_num, 1);

for step_id = 1:step_num
    vert_cell = cell(quality_poly_num, 1);
    tri_cell = cell(quality_poly_num, 1);
    vert_counter = 0;

    for poly_id = 1:quality_poly_num
        angle = angles(poly_id);
        rot_mat = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        vert_cell{poly_id} = sq_grid_cell{step_id} * rot_mat ...
            + shifts(poly_id, :);
        tri_cell{poly_id} = vert_counter + sq_tri_cell{step_id};
        vert_counter = vert_counter + size(vert_cell{poly_id}, 1);
    end

    quality_testfuncs_cell{step_id} = {vertcat(vert_cell{:}), ...
        vertcat(tri_cell{:}), false};
    quality_testfuncs_num_list(step_id) = vert_counter;
end

MT = MT2DQuad_MMOT(marg_cell, marg_weights, quality_poly_cell);

plot_quality_x_num = 3000;
plot_quality_y_num = 3000;
plot_quality_pt_x = linspace(quality_x_min, quality_x_max, ...
    plot_quality_x_num + 1);
plot_quality_pt_x = (plot_quality_pt_x(1:end - 1) ...
    + plot_quality_pt_x(2:end)) / 2;
plot_quality_pt_y = linspace(quality_y_min, quality_y_max, ...
    plot_quality_y_num + 1);
plot_quality_pt_y = (plot_quality_pt_y(1:end - 1) ...
    + plot_quality_pt_y(2:end)) / 2;
[plot_quality_grid_x, plot_quality_grid_y] = meshgrid( ...
    plot_quality_pt_x, plot_quality_pt_y);
plot_quality_pts = [plot_quality_grid_x(:), ...
    plot_quality_grid_y(:)];
plot_quality_inside = MT.checkIfInsideQualitySpace(plot_quality_pts);

quality_hist_x_num = 200;
quality_hist_y_num = 200;
quality_hist_edge_x = linspace(quality_x_min, quality_x_max, ...
    quality_hist_x_num + 1);
quality_hist_pt_x = (quality_hist_edge_x(1:end - 1) ...
    + quality_hist_edge_x(2:end)) / 2;
quality_hist_edge_y = linspace(quality_y_min, quality_y_max, ...
    quality_hist_y_num + 1);
quality_hist_pt_y = (quality_hist_edge_y(1:end - 1) ...
    + quality_hist_edge_y(2:end)) / 2;
[quality_hist_grid_x, quality_hist_grid_y] = meshgrid( ...
    quality_hist_pt_x, quality_hist_pt_y);
quality_hist_inside = MT.checkIfInsideQualitySpace( ...
    [quality_hist_grid_x(:), quality_hist_grid_y(:)]);

save(save_file_path, ...
    'marg_num', 'marg_weights', ...
    'marg_vert_cell', 'marg_tri_cell', 'marg_dens_cell', ...
    'mesh_sizes', 'testfuncs_cell', 'testfunc_num_mat', ...
    'quality_poly_cell', 'quality_testfuncs_cell', ...
    'quality_testfuncs_num_list', 'plot_pt_x_cell', 'plot_pt_y_cell', ...
    'plot_inside_cell', 'plot_pts_cell', ...
    'plot_quality_pt_x', 'plot_quality_pt_y', ...
    'plot_quality_inside', 'plot_quality_pts', ...
    'quality_hist_edge_x', 'quality_hist_edge_y', ...
    'quality_hist_pt_x', 'quality_hist_pt_y', 'quality_hist_inside');