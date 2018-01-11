function [center, Cn, PNR] = initComponents_residual_parallel(obj, K, save_avi, use_parallel, min_corr, min_pnr, seed_method)
%% initializing spatial/temporal components for the residual video
%% input:
%   K:  scalar, maximum number of neurons
%   frame_range: 1 X 2 vector indicating the starting and ending frames
%   save_avi: save the video of initialization procedure
%   use_parallel: boolean, do initialization in patch mode or not.
%       default(true); we recommend you to set it false only when you want to debug the code.

%% Output:
%   center: d*2 matrix, centers of all initialized neurons.
%   Cn:     correlation image
%   PNR:    peak to noise ratio
%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% process parameters

try
    % map data
    mat_data = obj.P.mat_data;
    
    % folders and files for saving the results
    log_folder = obj.P.log_folder;
    log_file = obj.P.log_file;
    log_data_file = obj.P.log_data;
    log_data = matfile(log_data_file, 'Writable', true); %#ok<NASGU>
    
    % dimension of data
    dims = mat_data.dims;
    d1 = dims(1);
    d2 = dims(2);
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
fprintf('\n---------------PICK NEURONS FROM THE RESIDUAL----------------\n');

if exist('min_corr', 'var') && ~isempty(min_corr)
    obj.options.min_corr = min_corr;
end
if exist('min_pnr', 'var') && ~isempty(min_pnr)
    obj.options.min_pnr = min_pnr;
end
if exist('seed_method', 'var') && ~isempty(seed_method)
    obj.options.seed_method = seed_method;
end
% frames to be loaded for initialization
frame_range = obj.frame_range;
T = diff(frame_range) + 1;

% maximum neuron number in each patch
if (~exist('K', 'var')) || (isempty(K))
    % if K is not specified, use a very large number as default
    K = round((d1*d2));
end

% exporting initialization procedures as a video
if ~exist('save_avi', 'var')||isempty(save_avi)
    save_avi = false; %don't save initialization procedure
elseif save_avi
    use_parallel = false;
end

% use parallel or not
if ~exist('use_parallel', 'var')||isempty(use_parallel)
    use_parallel = true; %don't save initialization procedure
end

% parameter for avoiding using boundaries pixels as seed pixels
options = obj.options;
if ~isfield(options, 'bd') || isempty(options.bd')
    options.bd = options.gSiz;   % boundary pixesl to be ignored during the process of detecting seed pixels
end
bd = options.bd;

if strcmpi(obj.options.seed_method, 'manual')
    use_parallel = false;
end
% no centering of the raw video
% options.center_psf = false;

%% prepare for the variables for computing the background.
bg_model = obj.options.background_model;
bg_ssub = obj.options.bg_ssub;
W = obj.W;
b0 = obj.b0;
b = obj.b;
f = obj.f;

%% the extracted neurons' signals
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

%% start initialization
% save the log infomation
k_options = obj.P.k_options +1;
eval(sprintf('log_data.options_%d=options; ', k_options));
obj.P.k_options = k_options;

flog = fopen(log_file, 'a');
fprintf(flog, '[%s]\b', get_minute());
fprintf(flog, 'Start picking neurons from the residual......\n\tThe collection of options are saved as intermediate_results.options_%d\n', k_options);

Ain = cell(nr_patch, nc_patch); % save spatial components of neurons in each patch
Cin = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, denoised trace
Sin = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, deconvolved trace
Cin_raw = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, raw trace
kernel_pars = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch
center = cell(nr_patch, nc_patch);     % save centers of all initialized neurons
Cn = zeros(d1, d2);
PNR = zeros(d1, d2);
default_kernel = obj.kernel;

results = cell(nr_patch*nc_patch, 1);
if use_parallel
    parfor mpatch=1:(nr_patch*nc_patch)
        % get the indices corresponding to the selected patch
        tmp_patch = patch_pos{mpatch};
        if strcmpi(bg_model, 'ring')
            % when it's ring model, log pixels surrounding the patch
            tmp_block = block_pos{mpatch};
        else
            tmp_block = tmp_patch;
        end
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        
        % use ind_patch to indicate pixels within the patch
        ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
        ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;
        
        % get the neural activity
        C_patch = C{mpatch};                % previous estimation of neural activity
        if isempty(C_patch)
            C_patch = 0;
            A_patch = 0;
        else
            A_patch = A{mpatch};
        end
        
        % boundaries pixels to be avoided for detecting seed pixels
        tmp_options = options;
        tmp_options.visible_off = true;
        tmp_options.bd = bd*([r==1, r==nr_patch, c==1, c==nc_patch]);
        
        % patch dimension
        tmp_options.d1 = diff(tmp_patch(1:2))+1;
        tmp_options.d2 = diff(tmp_patch(3:4))+1;
        
        % parameter for calcium indicators. This one may not be used if the
        % selected deconvolution algorithm doesn't need it
        tmp_options.kernel = default_kernel;
        
        % file names for saving avi file
        if save_avi
            tmp_save_avi = sprintf('%sinitialization_res_%d_%d_%d_%d.avi', log_folder, tmp_block(1), tmp_block(2), tmp_block(3), tmp_block(4));
        else
            tmp_save_avi = save_avi;
        end
        % get data
        if strcmpi(bg_model, 'ring')
            % including areas outside of the patch for recorving background
            % in the ring model
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        else
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, false);
        end
        [nr_block, nc_block, ~] = size(Ypatch);
        
        Ypatch = double(reshape(Ypatch, [], T)) - A_patch*C_patch;
        % get background
        if strcmpi(bg_model, 'ring')
            W_ring = W{mpatch};
            b0_ring = b0{mpatch};
            %             Ypatch = bsxfun(@minus, Ypatch(ind_patch,:)- W_ring*Ypatch, b0_ring-W_ring*mean(Ypatch, 2));
            if bg_ssub==1
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch,:))- W_ring*Ypatch, b0_ring-W_ring*mean(Ypatch, 2));
            else
                % get the dimension of the downsampled data
                [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
                % downsample data and reconstruct B^f
                temp = reshape(bsxfun(@minus, Ypatch, mean(Ypatch, 2)), nr_block, nc_block, []);
                temp = imresize(temp, 1./bg_ssub);
                Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
                Bf = imresize(Bf, [nr_block, nc_block]);
                Bf = reshape(Bf, [], T);
                
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch, :)) - Bf(ind_patch, :), b0_ring);
            end
        elseif strcmpi(bg_model, 'nmf')
            b_nmf = b{mpatch};
            f_nmf = f{mpatch};
            Ypatch = Ypatch- b_nmf*f_nmf;
        else
            b_svd = b{mpatch};
            f_svd = f{mpatch};
            b0_svd = b0{mpatch};
            Ypatch = Ypatch - bsxfun(@plus, b_svd*f_svd, b0_svd);
        end
        
        
        
        [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch, K, tmp_options, [], tmp_save_avi);
        
        % put everthing into one struct variable
        tmp_results.center = tmp_center;
        tmp_results.Cn = tmp_Cn;
        tmp_results.PNR = tmp_PNR;
        results{mpatch} = tmp_results;
        %     eval(sprintf('results_patch_%d=tmp_results;', mpatch));  %#ok<PFBFN>
        % display initialization progress
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
else
    for mpatch=1:(nr_patch*nc_patch)
        % get the indices corresponding to the selected patch
        tmp_patch = patch_pos{mpatch};
        if strcmpi(bg_model, 'ring')
            % when it's ring model, log pixels surrounding the patch
            tmp_block = block_pos{mpatch};
        else
            tmp_block = tmp_patch;
        end
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        
        % use ind_patch to indicate pixels within the patch
        ind_patch = false(diff(tmp_block(1:2))+1, diff(tmp_block(3:4))+1);
        ind_patch((tmp_patch(1):tmp_patch(2))-tmp_block(1)+1, (tmp_patch(3):tmp_patch(4))-tmp_block(3)+1) = true;
        
        % get the neural activity
        C_patch = C{mpatch};                % previous estimation of neural activity
        if isempty(C_patch)
            C_patch = 0;
            A_patch = 0;
        else
            A_patch = A{mpatch};
        end
        
        % boundaries pixels to be avoided for detecting seed pixels
        tmp_options = options;
        tmp_options.visible_off = true;
        tmp_options.bd = bd*([r==1, r==nr_patch, c==1, c==nc_patch]);
        
        % patch dimension
        tmp_options.d1 = diff(tmp_patch(1:2))+1;
        tmp_options.d2 = diff(tmp_patch(3:4))+1;
        
        % parameter for calcium indicators. This one may not be used if the
        % selected deconvolution algorithm doesn't need it
        tmp_options.kernel = default_kernel;
        
        % file names for saving avi file
        if save_avi
            tmp_save_avi = sprintf('%sinitialization_res_%d_%d_%d_%d.avi', log_folder, tmp_block(1), tmp_block(2), tmp_block(3), tmp_block(4));
        else
            tmp_save_avi = save_avi;
        end
        % get data
        if strcmpi(bg_model, 'ring')
            % including areas outside of the patch for recorving background
            % in the ring model
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        else
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, false);
        end
        [nr_block, nc_block, ~] = size(Ypatch);
        
        Ypatch = double(reshape(Ypatch, [], T)) - A_patch*C_patch;
        % get background
        if strcmpi(bg_model, 'ring')
            W_ring = W{mpatch};
            b0_ring = b0{mpatch};
            %             Ypatch = bsxfun(@minus, Ypatch(ind_patch,:)- W_ring*Ypatch, b0_ring-W_ring*mean(Ypatch, 2));
            if bg_ssub==1
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch,:))- W_ring*Ypatch, b0_ring-W_ring*mean(Ypatch, 2));
            else
                % get the dimension of the downsampled data
                [d1s, d2s] = size(imresize(zeros(nr_block, nc_block), 1/bg_ssub));
                % downsample data and reconstruct B^f
                temp = reshape(bsxfun(@minus, Ypatch, mean(Ypatch, 2)), nr_block, nc_block, []);
                temp = imresize(temp, 1./bg_ssub);
                Bf = reshape(W_ring*reshape(temp, [], T), d1s, d2s, T);
                Bf = imresize(Bf, [nr_block, nc_block]);
                Bf = reshape(Bf, [], T);
                
                Ypatch = bsxfun(@minus, double(Ypatch(ind_patch, :)) - Bf(ind_patch, :), b0_ring);
            end
            
        elseif strcmpi(bg_model, 'nmf')
            b_nmf = b{mpatch};
            f_nmf = f{mpatch};
            Ypatch = Ypatch- b_nmf*f_nmf;
        else
            b_svd = b{mpatch};
            f_svd = f{mpatch};
            b0_svd = b0{mpatch};
            Ypatch = Ypatch - bsxfun(@plus, b_svd*f_svd, b0_svd);
        end
        
        
        
        [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch, K, tmp_options, [], tmp_save_avi);
        
        % put everthing into one struct variable
        tmp_results.center = tmp_center;
        tmp_results.Cn = tmp_Cn;
        tmp_results.PNR = tmp_PNR;
        results{mpatch} = tmp_results;
        %     eval(sprintf('results_patch_%d=tmp_results;', mpatch));  %#ok<PFBFN>
        % display initialization progress
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
end

%% collect results
for mpatch=1:(nr_patch*nc_patch)
    % get the indices corresponding to the selected patch
    [mr, mc] = ind2sub([nr_patch, nc_patch], mpatch);
    tmp_patch = patch_pos{mr, mc};
    tmp_block = tmp_patch;
    r0 = tmp_block(1);
    r1 = tmp_block(2);
    c0 = tmp_block(3);
    c1 = tmp_block(4);
    %         ind_patch = true(r1-r0+1, c1-c0+1);
    %         ind_patch((tmp_patch(1):tmp_patch(2))-r0+1, (tmp_patch(3):tmp_patch(4))-c0+1) = false;
    %
    % unpack results
    tmp_results = results{mpatch};
    if isempty(tmp_results)
        continue;
    end
    
    tmp_Ain = tmp_results.Ain;
    tmp_ind = (sum(tmp_Ain, 1)>0);
    tmp_Ain = tmp_Ain(:, tmp_ind);
    tmp_Cin = tmp_results.Cin(tmp_ind,:);
    tmp_Cin_raw = tmp_results.Cin_raw(tmp_ind,:);
    tmp_center = tmp_results.center(tmp_ind,:) ;
    tmp_Cn = tmp_results.Cn;
    tmp_PNR = tmp_results.PNR;
    if options.deconv_flag
        tmp_Sin = tmp_results.Sin(tmp_ind,:);
        tmp_kernel_pars = tmp_results.kernel_pars(tmp_ind,:);
    end
    tmp_K = size(tmp_Ain, 2);   % number of neurons within the selected patch
    [tmp_d1, tmp_d2] = size(tmp_Cn);
    
    temp = zeros(d1, d2, tmp_K);  % spatial components of all neurons
    temp(r0:r1, c0:c1, :) = reshape(full(tmp_Ain), tmp_d1, tmp_d2, []);
    Ain{mr, mc} = reshape(temp, d1*d2, tmp_K);
    Cin{mr, mc} = tmp_Cin;      % temporal components of all neurons
    Cin_raw{mr, mc} = tmp_Cin_raw;
    if options.deconv_flag
        Sin{mr, mc} = tmp_Sin;
        kernel_pars{mr,mc} = tmp_kernel_pars;
    end
    center{mr, mc} = bsxfun(@plus, tmp_center, [r0-1, c0-1]); % centers
    
    Cn(r0:r1, c0:c1) = max(Cn(r0:r1, c0:c1), tmp_Cn);
    PNR(r0:r1, c0:c1) = max(PNR(r0:r1, c0:c1), tmp_PNR);
end
%% export the results
Ain = cell2mat(reshape(Ain, 1, []));
Cin = cell2mat(reshape(Cin, [], 1));
Cin_raw = cell2mat(reshape(Cin_raw, [], 1));
center = cell2mat(reshape(center, [], 1));
obj.A = [obj.A, sparse(Ain)];
obj.C = [obj.C; Cin];
obj.C_raw = [obj.C_raw; Cin_raw];
if options.deconv_flag
    Sin = cell2mat(reshape(Sin, [], 1));
    kernel_pars = cell2mat(reshape(kernel_pars, [], 1));
    obj.S = [obj.S; sparse(Sin)];
    obj.P.kernel_pars = [obj.P.kernel_pars; kernel_pars];
else
    obj.S = [obj.S; zeros(size(obj.C))];
end
K = size(Ain, 2);
k_ids = obj.P.k_ids;
obj.ids = [obj.ids, k_ids+(1:K)];
obj.tags = [obj.tags;zeros(K,1, 'like', uint16(0))];
obj.P.k_ids = K+k_ids;

%% save the results to log

fprintf(flog, '[%s]\b', get_minute());
fprintf(flog, 'Finished the initialization of neurons from the residual video.\n');
fprintf(flog, '\tIn total, %d neurons were detected. \n', size(Ain,2));

if obj.options.save_intermediate
    initialization_res.Ain = sparse(Ain);
    initialization_res.Cin = Cin;
    initialization_res.Cin_raw = Cin_raw;
    initialization_res.options = options;  %#ok<STRNU>
    
    tmp_str = get_date();
    tmp_str=strrep(tmp_str, '-', '_');
    eval(sprintf('log_data.initialization_res_%s = initialization_res;', tmp_str));
    
    fprintf(flog, '\tThe  results were saved as intermediate_results.initialization_res_%s\n\n', tmp_str);
end
fclose(flog);