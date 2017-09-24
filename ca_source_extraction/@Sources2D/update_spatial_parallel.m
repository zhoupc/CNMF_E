function update_spatial_parallel(obj, use_parallel)
%% update the the spatial components for all neurons
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% process parameters

try
    % map data
    mat_data = obj.P.mat_data;
    
    % folders and files for saving the results
    log_file =  obj.P.log_file;
    flog = fopen(log_file, 'a');
    log_data = matfile(obj.P.log_data, 'Writable', true); %#ok<NASGU>
    
    % dimension of data
    dims = mat_data.dims;
    d1 = dims(1);
    d2 = dims(2);
    T = dims(3);
    obj.options.d1 = d1;
    obj.options.d2 = d2;
    
    % parameters for patching information
    patch_pos = mat_data.patch_pos;
    block_pos = mat_data.block_pos;
    
    % number of patches
    [nr_patch, nc_patch] = size(patch_pos);
catch
    error('No data file selected');
end
fprintf('\n-----------------UPDATE SPATIAL---------------------------\n');
% frames to be loaded
frame_range = obj.frame_range;
T = diff(frame_range) + 1;

% threshold for detecting large residuals
thresh_outlier = obj.options.thresh_outlier;

% use parallel or not
if ~exist('use_parallel', 'var')||isempty(use_parallel)
    use_parallel = true; %don't save initialization procedure
end

% options
options = obj.options;
bg_model = options.background_model;

%% determine search location for each neuron
search_method = options.search_method;
if strcmpi(search_method, 'dilate')
    obj.options.se = [];
end
IND = sparse(logical(determine_search_location(obj.A, search_method, options)));

%% identify existing neurons within each patch
A = cell(nr_patch, nc_patch);
C = cell(nr_patch, nc_patch);
sn = cell(nr_patch, nc_patch);
ind_neurons = cell(nr_patch, nc_patch);
IND_search = cell(nr_patch, nc_patch);

for mpatch=1:(nr_patch*nc_patch)
    tmp_patch = patch_pos{mpatch};
    tmp_block = block_pos{mpatch};
    
    % find the neurons that are within the block
    if strcmpi(bg_model, 'ring')
        pause;
    else
        mask = false(d1, d2);
        mask(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4)) = true;
        ind = (reshape(double(mask(:)), 1, [])*IND>0);
        IND_search{mpatch} = IND(logical(mask), ind);
        A{mpatch}= obj.A(mask, ind);
        sn{mpatch} = obj.P.sn(mask);
        C{mpatch} = obj.C(ind, :);
        ind_neurons{mpatch} = find(ind);    % indices of the neurons within each patch
    end
    
end

%% prepare for the variables for computing the background.
if strcmpi(bg_model, 'ring')
    W = obj.W;
    b0 = obj.b0;
elseif strcmpi(bg_model, 'nmf')
    b = obj.b;
    f = obj.f;
else
    b = obj.b;
    f = obj.f;
    b0 = obj.b0;
end

%% start updating temporal components
A_new = A;
if use_parallel
    tmp_obj = Sources2D();
    tmp_obj.options = obj.options; 
    parfor mpatch=1:(nr_patch*nc_patch)
        % no neurons within the patch
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        tmp_patch = patch_pos{mpatch};     %[r0, r1, c0, c1], patch location
        A_patch = A{mpatch};
        if isempty(A_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        C_patch = C{mpatch};                % previous estimation of neural activity
        IND_patch = IND_search{mpatch};
        nr = diff(tmp_patch(1:2)) + 1;
        nc = diff(tmp_patch(3:4)) + 1;
        
        % get data
        if strcmpi(bg_model, 'ring')
            % including areas outside of the patch for recorving background
            % in the ring model
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        else
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, false);
        end
        
        % get background
        if strcmpi(bg_model, 'ring')
            pause;
        elseif strcmpi(bg_model, 'nmf')
            pause;
        else
            b_svd = b{mpatch};
            f_svd = f{mpatch};
            b0_svd = b0{mpatch};
            Ypatch = double(reshape(Ypatch, [], T)) - bsxfun(@plus, b_svd*f_svd, b0_svd);
        end
        
        % using HALS to update spatial components
        temp = HALS_spatial(Ypatch, A_patch, C_patch, IND_patch, 3);
        A_new{mpatch} = tmp_obj.post_process_spatial(reshape(full(temp), nr, nc, [])); %#ok<PFBNS>
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
else
    for mpatch=1:(nr_patch*nc_patch)
        % no neurons within the patch
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        tmp_patch = patch_pos{mpatch};     %[r0, r1, c0, c1], patch location
        A_patch = A{mpatch};
        if isempty(A_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        C_patch = C{mpatch};                % previous estimation of neural activity
        IND_patch = IND_search{mpatch};
        nr = diff(tmp_patch(1:2)) + 1;
        nc = diff(tmp_patch(3:4)) + 1;
        
        % get data
        if strcmpi(bg_model, 'ring')
            % including areas outside of the patch for recorving background
            % in the ring model
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        else
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, false);
        end
        
        % get background
        if strcmpi(bg_model, 'ring')
            pause;
        elseif strcmpi(bg_model, 'nmf')
            pause;
        else
            b_svd = b{mpatch};
            f_svd = f{mpatch};
            b0_svd = b0{mpatch};
            Ypatch = double(reshape(Ypatch, [], T)) - bsxfun(@plus, b_svd*f_svd, b0_svd);
        end
        
        % using HALS to update spatial components
        temp = HALS_spatial(Ypatch, A_patch, C_patch, IND_patch, 3);
        A_new{mpatch} = obj.post_process_spatial(reshape(full(temp), nr, nc, []));
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
end

%% collect results
K = size(obj.A, 2);
A_ = zeros(d1, d2, K);
fprintf('Collect results from all small patches...\n');
for mpatch=1:(nr_patch*nc_patch)
    A_patch = A_new{mpatch};
    tmp_pos = patch_pos{mpatch};
    nr = diff(tmp_pos(1:2))+1;
    nc = diff(tmp_pos(3:4))+1;
    ind_patch = ind_neurons{mpatch};
    for m=1:length(ind_patch)
        k = ind_patch(m);
        A_(tmp_pos(1):tmp_pos(2), tmp_pos(3):tmp_pos(4), k) = reshape(A_patch(:, m), nr, nc, 1);
    end
end
A_old = sparse(obj.A);
A_new = sparse(obj.reshape(A_, 1));

%% post-process results 
fprintf('Post-process spatial components of all neurons...\n');
obj.A = obj.post_process_spatial(A_new);
fprintf('Done!\n');
fprintf('Done!\n');

%% upadte b0
fprintf('Update the constant baselines for all pixels..\n');
C_mean = mean(obj.C, 2);
db = (A_old-A_new)*C_mean;
db = obj.reshape(db, 2);
b0 = obj.b0;
for mpatch = 1:(nr_patch*nc_patch)
    tmp_patch = patch_pos{mpatch};     %[r0, r1, c0, c1], patch location
    db_patch = db(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4));
    b0{mpatch} = b0{mpatch} + reshape(db_patch, [],1);
end
obj.b0 = b0;
fprintf('Done!\n');

%% save the results to log
spatial.A = obj.A;
temporal.b0 = obj.b0;
tmp_str = get_date();
tmp_str=strrep(tmp_str, '-', '_');
eval(sprintf('log_data.spatial_%s = spatial;', tmp_str));

fprintf(flog, '[%s]\b', get_minute());
fprintf(flog, 'Finished updating spatial components.\n');
fprintf(flog, '\tThe results were saved as intermediate_results.spatial_%s\n\n', tmp_str);
fclose(flog);