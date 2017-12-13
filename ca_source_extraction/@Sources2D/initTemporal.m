function initTemporal(obj, frame_range, use_parallel)
%% initializing temporal components with known spatial components and background information
%% input:
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% process parameters
try
    % map data
    mat_data = obj.P.mat_data;
    mat_file = mat_data.Properties.Source;
    
    % dimension of data
    dims = mat_data.dims;
    d1 = dims(1);
    d2 = dims(2);
    T = dims(3);
    obj.options.d1 = d1;
    obj.options.d2 = d2;
    % frames to be loaded for initialization
    if ~exist('frame_range', 'var')
        frame_range = obj.frame_range;
    end
    if isempty(frame_range)
        frame_range = [1, T];
    else
        frame_range(frame_range<1) = 1;
        frame_range(frame_range>T) = T;
    end
    T = diff(frame_range) + 1;
    obj.frame_range = frame_range;
    
    % folders and files for saving the results
    tmp_dir = sprintf('%s%sframes_%d_%d%s', fileparts(mat_file),filesep, frame_range(1), frame_range(2), filesep);
    if ~exist(tmp_dir, 'dir')
        mkdir(tmp_dir);
    end
    log_folder = [tmp_dir,  'LOGS_', get_date(), filesep];
    log_file = [log_folder, 'logs.txt'];
    log_data_file = [log_folder, 'intermediate_results.mat'];
    obj.P.log_folder = log_folder;
    obj.P.log_file = log_file;
    obj.P.log_data = log_data_file;
    log_data = matfile(log_data_file, 'Writable', true);
    
    % create a folder for new log
    mkdir(log_folder);
    
    % parameters for patching information
    patch_pos = mat_data.patch_pos;
    block_pos = mat_data.block_pos;
    
    % number of patches
    [nr_patch, nc_patch] = size(patch_pos);
catch
    error('No data file selected');
end

% use parallel or not
if ~exist('use_parallel', 'var')||isempty(use_parallel)
    use_parallel = true; %don't save initialization procedure
end

% options
options = obj.options;
bg_model = options.background_model;

%% create a folder for saving log information
% save the log infomation
log_data.options_0=options;
obj.P.k_options = 1;
obj.P.k_neurons = 0;
flog = fopen(log_file, 'w');
fprintf(flog, 'Data: %s\n\n', mat_file);
fprintf(flog, '--------%s--------\n', get_date());
fprintf(flog, '[%s]\b', get_minute());
fprintf(flog, 'Start running source extraction......\nThe collection of options are saved as intermediate_results.options_0\n\n');

fprintf(flog, '[%s]\b', get_minute());
fprintf(flog, 'Start initializing neurons from frame %d to frame %d\n\n', frame_range(1), frame_range(2));
fprintf(flog, 'The spatial components and the background are known already\n');

fprintf('\n----------------- INITIALIZE TEMPORAL COMPONENTS --------------------\n');


%% identify existing neurons within each patch
K = size(obj.A, 2);
A = cell(nr_patch, nc_patch);
C = cell(nr_patch, nc_patch);
sn = cell(nr_patch, nc_patch);
ind_neurons = cell(nr_patch, nc_patch);
Ymean = cell(nr_patch, nc_patch);

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
    ind_neurons{mpatch} = find(ind);    % indices of the neurons within each patch
end

%% prepare variables for computing the background.
bg_model = obj.options.background_model;
W = obj.W;
b0 = obj.b0;
b = obj.b;
f = obj.f;

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
        if strcmpi(bg_model, 'ring')
            tmp_block = block_pos{mpatch};
        else
            tmp_block = patch_pos{mpatch};
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
        temp = mean(Ypatch, 3);
        Ymean{mpatch} = temp((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1);
        Ypatch = reshape(Ypatch, [], T);
        
        if isempty(A_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        AA{mpatch}= sum(A_patch(ind_patch,:).^2, 1);
        
        % get background
        if strcmpi(bg_model, 'ring')
            W_ring = W{mpatch};
            Ypatch = double(Ypatch(ind_patch,:))-W_ring*double(Ypatch);
            A_patch = A_patch(ind_patch,:)-W_ring*A_patch;
            [~, C_patch] = HALS_temporal(Ypatch, A_patch, [], 10);
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch, C_patch, 2, deconv_options);
        elseif strcmpi(bg_model, 'nmf')
            b_nmf = b{mpatch};
            k = size(A_patch, 2);
            [~, tmp_C] = HALS_temporal(double(Ypatch), [A_patch, b_nmf], [], 10, []);
            C_patch = tmp_C(1:k, :);
            f_nmf = tmp_C((k+1):end, :);
            f{mpatch} = f_nmf;
            Ypatch = double(Ypatch)- b_nmf*f_nmf;
            [~, C_patch] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 10,[]);
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 2, deconv_options);
        else
            b_svd = b{mpatch};
            b0_svd = mean(Ypatch, 2);
            k = size(A_patch, 2);
            tmp_A = full([A_patch, b_svd]);
            tmp_C = (tmp_A'*tmp_A)\(tmp_A'*double(Ypatch)-tmp_A'*b0_svd*ones(1,T));
            C_patch = tmp_C(1:k, :);
            f_svd = tmp_C((k+1):end, :);
            f{mpatch} = f_svd;
            
            b0{mpatch} = b0_svd;
            Ypatch = double(Ypatch) - bsxfun(@plus, b_svd*f_svd, b0_svd);
            [~, C_patch] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 10,[]);
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 2, deconv_options);
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
        temp = mean(Ypatch, 3);
        Ymean{mpatch} = temp((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1);
        Ypatch = reshape(Ypatch, [], T);
        
        if isempty(A_patch)
            fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
            continue;
        end
        AA{mpatch}= sum(A_patch(ind_patch,:).^2, 1);
        
        % get background
        if strcmpi(bg_model, 'ring')
            W_ring = W{mpatch};
            Ypatch = double(Ypatch(ind_patch,:))-W_ring*double(Ypatch);
            A_patch = A_patch(ind_patch,:)-W_ring*A_patch;
            [~, C_patch] = HALS_temporal(double(Ypatch), A_patch, [], 10);
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch, C_patch, 2, deconv_options);
        elseif strcmpi(bg_model, 'nmf')
            b_nmf = b{mpatch};
            k = size(A_patch, 2);
            [~, tmp_C] = HALS_temporal(Ypatch, [A_patch, b_nmf], [], 10, []);
            C_patch = tmp_C(1:k, :);
            f_nmf = tmp_C((k+1):end, :);
            f{mpatch} = f_nmf;
            Ypatch = double(Ypatch)- b_nmf*f_nmf;
            [~, C_patch] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 10,[]);
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 2, deconv_options);
        else
            b_svd = b{mpatch};
            b0_svd = mean(Ypatch, 2);
            k = size(A_patch, 2);
            tmp_A = [A_patch, b_svd];
            tmp_C = (tmp_A'*tmp_A)\(tmp_A'*double(Ypatch)-tmp_A'*b0_svd*ones(1,T));
            C_patch = tmp_C(1:k, :);
            f_svd = tmp_C((k+1):end, :);
            f{mpatch} = f_svd;
            
            b0{mpatch} = b0_svd;
            Ypatch = double(Ypatch) - bsxfun(@plus, b_svd*f_svd, b0_svd);
            [~, C_patch] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 10,[]);
            [~, C_raw_new{mpatch}] = HALS_temporal(Ypatch, A_patch(ind_patch,:), C_patch, 2, deconv_options);
        end
        
        
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
end

obj.b0 = b0;
obj.b=b;
obj.f=f;
obj.P.Ymean = Ymean;

%% collect results
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
obj.deconvTemporal();
fprintf('Done!\n');

%% save the results to log
fprintf(flog, '[%s]\b', get_minute());
if obj.options.save_intermediate
    initialization.neuron = obj.obj2struct();
    log_data.initialization = initialization;
    fprintf(flog, '\tThe initialization results were saved as intermediate_results.initialization\n\n');
end
fprintf(flog, 'Finished the initialization procedure.\n');
fclose(flog);