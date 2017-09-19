function update_temporal_parallel(obj, use_parallel)
%% update the the temporal components for all neurons
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
fprintf('\n-----------------UPDATE TEMPORAL---------------------------\n');

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

%% identify existing neurons within each patch
A = cell(nr_patch, nc_patch);
C = cell(nr_patch, nc_patch);
sn = cell(nr_patch, nc_patch);
ind_neurons = cell(nr_patch, nc_patch);

AA = cell(nr_patch, nc_patch);   % save the ai^T*ai for each neuron

for mpatch=1:(nr_patch*nc_patch)
    if strcmpi(bg_model, 'ring')
        tmp_block = block_pos{mpatch};
    else
        tmp_block = patch_pos{mpatch};
    end
    % find the neurons that are within the block
    mask = zeros(d1, d2);
    mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
    ind = (reshape(mask(:), 1, [])* obj.A>0);
    A{mpatch}= obj.A(logical(mask), ind);
    AA{mpatch}= sum(A{mpatch}.^2, 1);
    sn{mpatch} = obj.P.sn(logical(mask));
    C{mpatch} = obj.C(ind, :);
    ind_neurons{mpatch} = find(ind);    % indices of the neurons within each patch
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
C_raw_new = C;
if obj.options.deconv_flag
    deconv_options = obj.options.deconv_options;
else
    deconv_options = [];
end
if use_parallel
    parfor mpatch=1:(nr_patch*nc_patch)
        % no neurons within the patch
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        tmp_patch = patch_pos{mpatch};     %[r0, r1, c0, c1], patch location
        C_patch = C{mpatch};                % previous estimation of neural activity
        if isempty(C_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        A_patch = A{mpatch};
        
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
        
        % using HALS to update temporal components
        [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch, C_patch, 2, deconv_options);
        
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
else
    for mpatch=1:(nr_patch*nc_patch)
        % no neurons within the patch
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        tmp_patch = patch_pos{mpatch};     %[r0, r1, c0, c1], patch location
        C_patch = C{mpatch};                % previous estimation of neural activity
        if isempty(C_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        A_patch = A{mpatch};
        
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
        
        % using HALS to update temporal components
        [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch, C_patch, 2, deconv_options);
        
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
end
%% collect results
K = size(obj.C, 1);
C_new = zeros(K, T);
aa = zeros(K, 1);
fprintf('Collect results from all small patches...\n');
for mpatch=1:(nr_patch*nc_patch)
    C_raw_patch = C_raw_new{mpatch};
    ind_patch = ind_neurons{mpatch};
    aa_patch = AA{mpatch};
    for m=1:length(ind_patch)
        k = ind_patch(m);
        C_new(k, :) = C_new(k, :) + C_raw_patch(m,:) * aa_patch(m);
        aa(k) = aa(k) + aa_patch(m);
    end
end
aa(aa==0) = 1;
obj.C_raw = bsxfun(@times, C_new, 1./aa);
fprintf('Deconvolve and denoise all temporal traces again...\n');
C_new = obj.deconvTemporal();
fprintf('Done!\n');

%% upadte b0
fprintf('Update the constant baselines for all pixels..\n');
C_old = obj.C;
db = obj.A*(mean(C_old, 2)-mean(C_new,2));
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
temporal.C_raw = obj.C_raw;
temporal.C = obj.C;
temporal.S = obj.S;
temporal.P.kernel_pars = obj.P.kernel_pars;
temporal.b0 = obj.b0;
tmp_str = get_date();
tmp_str=strrep(tmp_str, '-', '_');
eval(sprintf('log_data.temporal_%s = temporal;', tmp_str));

fprintf(flog, '%s\b', get_minute());
fprintf(flog, 'Finished updating temporal components.\n');
fprintf(flog, 'The results were saved as intermediate_results.temporal_%s\n\n', tmp_str);
fclose(flog);