function update_temporal_parallel(obj, use_parallel, use_c_hat)
%% update the the temporal components for all neurons
% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.
%   use_c_hat: use the previous estimation of C

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

% use previous estimation of C
if ~exist('use_c_hat', 'var') || isempty(use_c_hat)
    use_c_hat = true;
end
% options
options = obj.options;
bg_model = options.background_model;
maxIter = options.maxIter;

%% identify existing neurons within each patch
A = cell(nr_patch, nc_patch);
C = cell(nr_patch, nc_patch);
% if strcmpi(bg_model, 'ring')
    A_prev = A;
    C_prev = C;
% end
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
    sn{mpatch} = obj.P.sn(logical(mask));
    C{mpatch} = obj.C(ind, :);
    ind_neurons{mpatch} = find(ind);    % indices of the neurons within each patch
    
    if strcmpi(bg_model, 'ring')
        ind = find(reshape(mask(:)==1, 1, [])* full(obj.A_prev)>0);
        A_prev{mpatch}= obj.A_prev((mask>0), ind);
        C_prev{mpatch} = obj.C_prev(ind, :);
    end
end
%% prepare for the variables for computing the background.
bg_model = obj.options.background_model;
bg_ssub = obj.options.bg_ssub;
W = obj.W;
b0 = obj.b0;
b = obj.b;
f = obj.f;

%% start updating temporal components
C_raw_new = C;
deconv_flag = obj.options.deconv_flag;
if options.deconv_flag
    deconv_options = obj.options.deconv_options;
else
    deconv_options = [];
end
if use_parallel
    parfor mpatch=1:(nr_patch*nc_patch)
        % no neurons within the patch
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        tmp_patch = patch_pos{mpatch};     %[r0, r1, c0, c1], patch location
        if strcmpi(bg_model, 'ring')
            tmp_block = block_pos{mpatch};
        else
            tmp_block = patch_pos{mpatch};
        end
        C_patch = C{mpatch};                % previous estimation of neural activity
        
        if isempty(C_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        A_patch = A{mpatch};
        
        % use ind_patch to indicate pixels within the patch
        ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
        ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;
        
        % get data
        if strcmpi(bg_model, 'ring')
            % including areas outside of the patch for recorving background
            % in the ring model
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        else
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, false);
        end
        [nr_block, nc_block, ~] = size(Ypatch);
        
        % get background
        if strcmpi(bg_model, 'ring')
            A_patch_prev = A_prev{mpatch};
            C_patch_prev = C_prev{mpatch};
            W_ring = W{mpatch};
            b0_ring = b0{mpatch};
            Ypatch = reshape(Ypatch, [], T);
            tmp_Y = double(Ypatch)-A_patch_prev*C_patch_prev;
            if bg_ssub==1
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch,:))- W_ring*tmp_Y, b0_ring-W_ring*mean(tmp_Y, 2));
            else
                % get the dimension of the downsampled data
                [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
                % downsample data and reconstruct B^f
                temp = reshape(bsxfun(@minus, tmp_Y, mean(tmp_Y, 2)), nr_block, nc_block, []);
                temp = imresize(temp, 1./bg_ssub);
                Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
                Bf = imresize(Bf, [nr_block, nc_block]);
                Bf = reshape(Bf, [], T);
                
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch, :)) - Bf(ind_patch, :), b0_ring);
            end
        elseif strcmpi(bg_model, 'nmf')
            b_nmf = b{mpatch};
            f_nmf = f{mpatch};
            Ypatch = double(reshape(Ypatch, [], T))- b_nmf*f_nmf;
        else
            b_svd = b{mpatch};
            f_svd = f{mpatch};
            b0_svd = b0{mpatch};
            Ypatch = double(reshape(Ypatch, [], T)) - bsxfun(@plus, b_svd*f_svd, b0_svd);
        end
        
        % using HALS to update temporal components
        if ~use_c_hat
            [AA{mpatch}, C_raw_new{mpatch}] = fast_temporal(Ypatch, A_patch(ind_patch,:));
        else
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, maxIter, deconv_options);
            AA{mpatch}= sum(A_patch(ind_patch,:).^2, 1);
        end
        
        
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
else
    for mpatch=1:(nr_patch*nc_patch)
        % no neurons within the patch
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        tmp_patch = patch_pos{mpatch};     %[r0, r1, c0, c1], patch location
        if strcmpi(bg_model, 'ring')
            tmp_block = block_pos{mpatch};
        else
            tmp_block = patch_pos{mpatch};
        end
        
        C_patch = C{mpatch};                % previous estimation of neural activity
        
        if isempty(C_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        A_patch = A{mpatch};
        
        % use ind_patch to indicate pixels within the patch
        ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
        ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;
        
        % get data
        if strcmpi(bg_model, 'ring')
            % including areas outside of the patch for recorving background
            % in the ring model
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        else
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, false);
        end
        [nr_block, nc_block, ~] = size(Ypatch);
        
        % get background
        if strcmpi(bg_model, 'ring')
            A_patch_prev = A_prev{mpatch};
            C_patch_prev = C_prev{mpatch};
            W_ring = W{mpatch};
            b0_ring = b0{mpatch};
            Ypatch = reshape(Ypatch, [], T);
            tmp_Y = double(Ypatch)-A_patch_prev*C_patch_prev;
            if bg_ssub==1
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch,:))- W_ring*tmp_Y, b0_ring-W_ring*mean(tmp_Y, 2));
            else
                % get the dimension of the downsampled data
                [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
                % downsample data and reconstruct B^f
                temp = reshape(bsxfun(@minus, tmp_Y, mean(tmp_Y, 2)), nr_block, nc_block, []);
                temp = imresize(temp, 1./bg_ssub);
                Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
                Bf = imresize(Bf, [nr_block, nc_block]);
                Bf = reshape(Bf, [], T);
                
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch, :)) - Bf(ind_patch, :), b0_ring);
            end
        elseif strcmpi(bg_model, 'nmf')
            b_nmf = b{mpatch};
            f_nmf = f{mpatch};
            Ypatch = double(reshape(Ypatch, [], T))- b_nmf*f_nmf;
        else
            b_svd = b{mpatch};
            f_svd = f{mpatch};
            b0_svd = b0{mpatch};
            Ypatch = double(reshape(Ypatch, [], T)) - bsxfun(@plus, b_svd*f_svd, b0_svd);
        end
        
        % using HALS to update temporal components
        if ~use_c_hat
            [AA{mpatch}, C_raw_new{mpatch}] = fast_temporal(Ypatch, A_patch(ind_patch,:));
        else
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, maxIter, deconv_options);
            AA{mpatch}= sum(A_patch(ind_patch,:).^2, 1);
        end
        
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
if deconv_flag
    obj.C = obj.deconvTemporal();
else
    obj.C_raw = bsxfun(@minus, obj.C_raw, min(obj.C_raw,[],2)); 
    obj.C = obj.C_raw; 
end
fprintf('Done!\n');

%% upadte b0
if strcmpi(bg_model, 'ring')
    fprintf('Update the constant baselines for all pixels..\n');
    obj.b0_new = cell2mat(obj.P.Ymean)-obj.reshape(obj.A*mean(obj.C,2), 2); -obj.reconstruct_b0();
    fprintf('Done!\n');
end

%% save the results to log
fprintf(flog, '[%s]\b', get_minute());
fprintf(flog, 'Finished updating temporal components.\n');
if obj.options.save_intermediate
    temporal.C_raw = obj.C_raw;
    temporal.ids = obj.ids;
    temporal.C = obj.C;
    temporal.S = obj.S;
    temporal.P.kernel_pars = obj.P.kernel_pars;
    temporal.b0_new = obj.b0_new;
    tmp_str = get_date();
    tmp_str=strrep(tmp_str, '-', '_');
    eval(sprintf('log_data.temporal_%s = temporal;', tmp_str));
    fprintf(flog, '\tThe results were saved as intermediate_results.temporal_%s\n\n', tmp_str);
end
fclose(flog);

function [aa, C_raw] = fast_temporal(Y, A)
%% estimate temporal components using the mean fluorescence
% input:
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components
% output:
%   aa: K*1, energy of each Ai
%   C_raw: K*T, the temporal components before thresholded or being
%   denoised.
% Author: Pengcheng Zhou, Columbia University, 2017
% zhoupc1988@gmail.com

%% options

%% initialization
% roughly initialize C
tmpA = bsxfun(@times, A, 1./max(A, [], 1));
ind_max = bsxfun(@ge, tmpA, 0.5); %max(tmpA, [], 2));
tmp_A = A.*double(ind_max);
aa = sum(tmp_A.^2, 1);
ind = (aa==0);
aa(ind) = inf;
C_raw = bsxfun(@times, tmp_A'*Y, 1./aa');
aa(ind) = 0;
% if any(ind)  % explain the residual using weak neurons
%     tmp_A = A(:, ind);
%     C_raw(ind, :) = (tmp_A'*tmp_A)\(tmp_A'*Y - tmp_A'*A*C_raw);
%     aa(ind) = sum(tmp_A.^2, 1);
% end













