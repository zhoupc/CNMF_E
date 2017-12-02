function update_background_parallel(obj, use_parallel)
%% update the background related variables in CNMF framework
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
fprintf('\n-----------------UPDATE BACKGROUND---------------------------\n');

% frames to be loaded for initialization
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
nb = options.nb;
bg_ssub = options.bg_ssub;
bg_model = options.background_model;
with_projection = options.bg_acceleration;

% previous estimation
A = cell(nr_patch, nc_patch);
C = cell(nr_patch, nc_patch);
sn = cell(nr_patch, nc_patch);
W = obj.W;
b0 = obj.b0;
b = obj.b;
f = obj.f;
RSS = cell(nr_patch, nc_patch);

%% check whether the bg_ssub was changed
if strcmpi(bg_model, 'ring')
    tmp_block = block_pos{1};    % block position
    nr_block = diff(tmp_block(1:2))+1;
    nc_block = diff(tmp_block(3:4))+1;
    [~, temp] = size(W{1});
    [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
    if temp~=d1s*d2s
        rr = ceil(obj.options.ring_radius/bg_ssub);    % radius of the ring
        [r_shift, c_shift] = get_nhood(rr, obj.options.num_neighbors);    % shifts used for acquiring the neighboring pixels on the ring
        parfor mpatch=1:(nr_patch*nc_patch)
            tmp_patch = patch_pos{mpatch};    % patch position
            tmp_block = block_pos{mpatch};    % block position
            nr = diff(tmp_patch(1:2)) + 1;
            nc = diff(tmp_patch(3:4)) + 1;
            nr_block = diff(tmp_block(1:2))+1;
            nc_block = diff(tmp_block(3:4))+1;
            b0{mpatch} = zeros(nr*nc, 1);
            
            if bg_ssub==1
                [csub, rsub] = meshgrid(tmp_patch(3):tmp_patch(4), tmp_patch(1):tmp_patch(2));
                csub = reshape(csub, [], 1);
                rsub = reshape(rsub, [], 1);
                ii = repmat((1:numel(csub))', [1, length(r_shift)]);
                csub = bsxfun(@plus, csub, c_shift);
                rsub = bsxfun(@plus, rsub, r_shift);
                ind = and(and(csub>=1, csub<=d2), and(rsub>=1, rsub<=d1));
                jj = (csub-tmp_block(3)) * (diff(tmp_block(1:2))+1) + (rsub-tmp_block(1)+1);
                
                temp = sparse(ii(ind), jj(ind), 1, nr*nc, nr_block*nc_block);
                W{mpatch} = bsxfun(@times, temp, 1./sum(temp, 2));
            else
                d1s = ceil(nr_block/bg_ssub);
                d2s = ceil(nc_block/bg_ssub);
                
                [csub, rsub] = meshgrid(1:d2s, 1:d1s);
                csub = reshape(csub, [], 1);
                rsub = reshape(rsub, [], 1);
                ii = repmat((1:numel(csub))', [1, length(r_shift)]);
                csub = bsxfun(@plus, csub, c_shift);
                rsub = bsxfun(@plus, rsub, r_shift);
                jj = (csub-1) * d1s + rsub;
                % remove neighbors that are out of boundary
                ind = and(and(csub>=1, csub<=d2s), and(rsub>=1, rsub<=d1s));
                temp = sparse(ii(ind), jj(ind), 1, d1s*d2s, d1s*d2s);
                W{mpatch} = bsxfun(@times, temp, 1./sum(temp, 2));
            end
        end
    end
end

%% start updating the background
for mpatch=1:(nr_patch*nc_patch)
    tmp_block = block_pos{mpatch};
    
    % find the neurons that are within the block
    mask = zeros(d1, d2);
    mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1;
    
    ind = (reshape(mask(:), 1, [])* obj.A>0);
    A{mpatch}= obj.A(logical(mask), ind);
    C{mpatch} = obj.C(ind, :);
    temp = obj.P.sn(logical(mask));
    if bg_ssub==1
        sn{mpatch} = temp;
    else
        nr_block = diff(tmp_block(1:2))+1;
        nc_block = diff(tmp_block(3:4))+1;
        sn{mpatch} = imresize(reshape(temp, nr_block, nc_block), 1/bg_ssub, 'nearest')*bg_ssub;
    end
end

% check whether this is the first run of updating background components
if strcmpi(bg_model, 'ring')
    flag_first = (length(unique(W{1}(1, :)))==2);
else
    flag_first = (mean2(b{1})==0);
end

if use_parallel
    % load data before running parfor
    %     data_patch = cell(nr_patch, nc_patch);
    %     flag_ignore = cell(nr_patch, nc_patch);
    %     for mpatch=1:(nr_patch*nc_patch)
    %         tmp_patch = patch_pos{mpatch};
    %         A_block = A{mpatch};
    %
    %         % stop updating B because A&C doesn't change in this area
    %         if isempty(A_block) && (~flag_first)
    %             flag_ignore{mpatch} = true;
    %             [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
    %
    %             % keep the current results. this step looks rediculous, but it
    %             % is needed for some computer/matlab. very weird.
    %             W{mpatch} = W{mpatch};
    %             b0{mpatch} = b0{mpatch};
    %             b{mpatch} = b{mpatch};
    %             f{mpatch} = f{mpatch};
    %             fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    %             continue;
    %         else
    %             flag_ignore{mpatch} = false;
    %             % pull data
    %             data_patch{mpatch} = get_patch_data(mat_data, tmp_patch, frame_range, true);
    %         end
    %     end
    % do the actual computation
    parfor mpatch=1:(nr_patch*nc_patch)
        %         if flag_ignore{mpatch}
        %             continue;
        %         end
        tmp_patch = patch_pos{mpatch};
        tmp_block = block_pos{mpatch};
        
        A_block = A{mpatch};
        sn_block = sn{mpatch};
        C_block = C{mpatch};
        
        % stop updating B because A&C doesn't change in this area
        if isempty(A_block) && (~flag_first)
            [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
            
            % keep the current results. this step looks rediculous, but it
            % is needed for some computer/matlab. very weird.
            W{mpatch} = W{mpatch};
            b0{mpatch} = b0{mpatch};
            b{mpatch} = b{mpatch};
            f{mpatch} = f{mpatch};
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        
        % use ind_patch to indicate pixels within the patch and only
        % update (W, b0) corresponding to these pixels
        ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
        ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;
        
        % pull data
        %         Ypatch = data_patch{mpatch};
        Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        [nr_block, nc_block, T_block] = size(Ypatch);
        if strcmpi(bg_model, 'ring')
            % get the previous estimation
            W_old = W{mpatch};
            Ypatch = reshape(Ypatch, [], T);
            
            % run regression to get A, C, and W, b0
            if bg_ssub==1
                sn_patch = sn_block(ind_patch);
                [W{mpatch}, b0{mpatch}] = fit_ring_model(Ypatch, A_block, C_block, W_old, thresh_outlier, sn_patch, ind_patch, with_projection);
            else
                % downsapmle data first
                temp = reshape(double(Ypatch)-A_block*C_block, nr_block, nc_block, T_block);
                tmp_b0 = mean(temp, 3);
                b0{mpatch} = tmp_b0(ind_patch);
                Ypatch = imresize(temp, 1./bg_ssub, 'nearest');
                Ypatch = reshape(Ypatch, [], T_block);
                
                [W{mpatch}, ~] = fit_ring_model(Ypatch, [], [], W_old, thresh_outlier, sn_block(:), [],  with_projection);
                %                 tmp_b0 = imresize(reshape(tmp_b0, size(sn_block)), [nr_block, nc_block]);
                %                 b0{mpatch} = tmp_b0(ind_patch(:));
            end
        elseif strcmpi(bg_model, 'nmf')
            b_old = b{mpatch};
            f_old = f{mpatch};
            Ypatch = reshape(Ypatch, [], T);
            sn_patch = sn_block(ind_patch);
            [b{mpatch}, f{mpatch}] = fit_nmf_model(Ypatch, nb, A_block, C_block, b_old, f_old, thresh_outlier, sn_patch, ind_patch);
        else
            b_old = b{mpatch};
            f_old = f{mpatch};
            Ypatch = reshape(Ypatch, [], T);
            sn_patch = sn_block(ind_patch);
            [b{mpatch}, f{mpatch}, b0{mpatch}] = fit_svd_model(Ypatch, nb, A_block, C_block, b_old, f_old, thresh_outlier, sn_patch, ind_patch);
        end
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
        
    end
else
    for mpatch=1:(nr_patch*nc_patch)
        tmp_patch = patch_pos{mpatch};
        tmp_block = block_pos{mpatch};
        
        A_block = A{mpatch};
        sn_block = sn{mpatch};
        C_block = C{mpatch};
        
        % stop the updating B because A&C doesn't change in this area
        if isempty(A_block) && (~flag_first)
            [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        
        % use ind_patch to indicate pixels within the patch and only
        % update (W, b0) corresponding to these pixels
        ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
        ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;
        
        % pull data
        Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        [nr_block, nc_block, T_block] = size(Ypatch);
        if strcmpi(bg_model, 'ring')
            % get the previous estimation
            W_old = W{mpatch};
            Ypatch = reshape(Ypatch, [], T_block);
            
            % run regression to get A, C, and W, b0
            if bg_ssub==1
                sn_patch = sn_block(ind_patch);
                [W{mpatch}, b0{mpatch}] = fit_ring_model(Ypatch, A_block, C_block, W_old, thresh_outlier, sn_patch, ind_patch, with_projection);
            else
                % downsapmle data first
                temp = reshape(double(Ypatch)-A_block*C_block, nr_block, nc_block, T_block);
                tmp_b0 = mean(temp, 3);
                b0{mpatch} = tmp_b0(ind_patch);
                Ypatch = imresize(temp, 1./bg_ssub, 'nearest');
                Ypatch = reshape(Ypatch, [], T_block);
                
                [W{mpatch}, ~] = fit_ring_model(Ypatch, [], [], W_old, thresh_outlier, sn_block(:), [],  with_projection);
                %                 tmp_b0 = imresize(reshape(tmp_b0, size(sn_block)), [nr_block, nc_block]);
                %                 b0{mpatch} = tmp_b0(ind_patch(:));
            end
        elseif strcmpi(bg_model, 'nmf')
            b_old = b{mpatch};
            f_old = f{mpatch};
            Ypatch = reshape(Ypatch, [], T);
            sn_patch = sn_block(ind_patch);
            [b{mpatch}, f{mpatch}] = fit_nmf_model(Ypatch, nb, A_block, C_block, b_old, f_old, thresh_outlier, sn_patch, ind_patch);
        else
            b_old = b{mpatch};
            f_old = f{mpatch};
            Ypatch = reshape(Ypatch, [], T);
            sn_patch = sn_block(ind_patch);
            [b{mpatch}, f{mpatch}, b0{mpatch}] = fit_svd_model(Ypatch, nb, A_block, C_block, b_old, f_old, thresh_outlier, sn_patch, ind_patch);
        end
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
        
    end
end
obj.b = b;
obj.f = f;
obj.b0 = b0;
obj.W = W;
obj.b0_new = obj.reconstruct_b0();
obj.A_prev = obj.A;
obj.C_prev = obj.C;

%% save the results to log
fprintf('Finished updating background using %s model.\n', bg_model);

fprintf(flog, '[%s]\b', get_minute());
fprintf(flog, 'Finished updating background using %s model.\n', bg_model);
if obj.options.save_intermediate
    bg.b = obj.b;
    bg.f = obj.f;
    bg.b0 = obj.b0;
    bg.W = obj.W;
    tmp_str = get_date();
    tmp_str=strrep(tmp_str, '-', '_');
    eval(sprintf('log_data.bg_%s = bg;', tmp_str));
    fprintf(flog, '\tThe results were saved as intermediate_results.bg_%s\n\n', tmp_str);
end
fclose(flog);