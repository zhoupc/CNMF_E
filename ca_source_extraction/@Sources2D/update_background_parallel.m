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
    mat_file = mat_data.Properties.Source;
    tmp_dir = fileparts(mat_file);
    
    % folders and files for saving the results
    log_folder = [tmp_dir, filesep, 'LOGS_', get_date(), filesep];
    log_file = [log_folder, 'logs.txt'];
    log_data_file = [log_folder, 'intermediate_results.mat'];
    obj.P.log_folder = log_folder;
    obj.P.log_file = log_file;
    obj.P.log_data = log_data_file;
    log_data = matfile(log_data_file, 'Writable', true);
    
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

% frames to be loaded for initialization
frame_range = obj.frame_range;
T = diff(frame_range) + 1;

% threshold for detecting large residuals
thresh_outier = obj.options.thresh_outlier; 

% use parallel or not 
if ~exist('use_parallel', 'var')||isempty(use_parallel)
    use_parallel = true; %don't save initialization procedure
end

% options 
options = obj.options;

%% start updating the background 
bg_model = obj.options.background_model; 

if strcmpi(bg_model, 'ring')
    if use_parallel 
        % update each patch 
        for mpatch=1:(nr_patch*nc_patch)
            tmp_patch = patch_pos{mpatch}; 
            tmp_block = block_pos{mpatch};
            
            % extract the old (W, b)
            W_old = obj.W{mpatch}; 
            b0_old = obj.b0{mpatch};
            
            % find the neurons that are within the block 
            mask = zeros(d1, d2); 
            mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4)) = 1; 
            ind = (reshape(mask(:), 1, [])*obj.A>0); 
            A_patch = obj.A(:, ind); 
            C_patch = obj.C(:, ind); 
            
            % use ind_patch to indicate pixels within the patch and only
            % update (W, b0) corresponding to these pixels
            mask(tmp_patch(1):tmp_patch(2), tmp_patch(3):tmp_patch(4)) = 0; 
            ind_patch = logical(1-mask(tmp_block(1):tmp_block(2), tmp_block(3):tmp_block(4))); 
            
            % pull data 
            Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
            
            % run regression to get A, C, and W, b0
            [obj.W{match}, obj.b0{match}, ~] = fit_ring_model(Ypatch, A_patch, C_patch, W_old, b0_old, thresh_outlier, ind_patch);
            
            
            
        end
%     else 
%         for ...
%         end 
    end 
    
    
end 
%% model used for updating the background 
bg_model = obj.options.background_model;
if strcmpi(bg_model, 'ring')   
    rr = obj.options.ring_radius;    % radius of the ring 
    [r_shift, c_shift] = get_nhood(rr);    % shifts used for acquiring the neighboring pixels on the ring 
    W = cell(nr_patch, nc_patch);    % matrix for saving the weight matrix within each block 
    b0 = cell(nr_patch, nc_patch);   % constant baselines for all pixels 
    
    parfor mpatch=1:(nr_patch*nc_patch)
        tmp_patch = patch_pos{mpatch};    % patch position 
        tmp_block = block_pos{mpatch};    % block position 
        nr = diff(tmp_patch(1:2)) + 1; 
        nc = diff(tmp_patch(3:4)) + 1; 
        nr_block = diff(tmp_block(1:2))+1; 
        nc_block = diff(tmp_block(3:4))+1; 
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
        b0{mpatch} = zeros(nr*nc, 1); 
    end
    obj.W = W; 
    obj.b0 = b0; 
    clear W b0; 
elseif strcmpi(bg_model, 'nmf')
    b = cell(nr_patch, nc_patch);
    f = cell(nr_patch, nc_patch);
    nb = obj.options.nb;
    for mpatch=1:(nr_patch*nc_patch)
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);   % patch ind
        tmp_patch = patch_pos{r, c};    % patch position
        nr = diff(tmp_patch(1:2)) + 1;
        nc = diff(tmp_patch(3:4)) + 1;
        b{r, c} = zeros(nr*nc, nb);
        f{r, c} = zeros(nb, T);
    end
    obj.b = b; 
    obj.f = f; 
else
    %default, SVD model
    b = cell(nr_patch, nc_patch);
    f = cell(nr_patch, nc_patch);
    b0 = cell(nr_patch, nc_patch); 
    nb = obj.options.nb;
    for mpatch=1:(nr_patch*nc_patch)
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);   % patch ind
        tmp_patch = patch_pos{r, c};    % patch position
        nr = diff(tmp_patch(1:2)) + 1;
        nc = diff(tmp_patch(3:4)) + 1;
        b{r, c} = zeros(nr*nc, nb);
        f{r, c} = zeros(nb, T);
        b0{r, c} = zeros(nr*nc, 1); 
    end
    obj.b = b; 
    obj.f = f; 
    obj.b0 = b0; 
end

%% start initialization
% save the log infomation
log_data.options_0=options;
obj.P.k_options = 1;
obj.P.k_neurons = 0;
flog = fopen(log_file, 'w');
fprintf(flog, 'Data: %s\n\n', mat_file);

fprintf(flog, '%s\b', get_minute());
fprintf(flog, 'Start running source extraction......\nThe collection of options are saved as intermediate_results.options_0\n\n');

fprintf(flog, '%s\b', get_minute());
fprintf(flog, 'Start initializing neurons from frame %d to frame %d\n\n', frame_range(1), frame_range(2));

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
        tmp_block = block_pos{mpatch}; 
        
        % boundaries pixels to be avoided for detecting seed pixels
        tmp_options = options;
        tmp_options.visible_off = true;
        tmp_options.bd = max([(tmp_patch-tmp_block).*[1, -1, 1, -1]; bd, bd, bd, bd], [], 1);
        
        % patch dimension
        tmp_options.d1 = diff(tmp_block(1:2))+1;
        tmp_options.d2 = diff(tmp_block(3:4))+1;
        
        % parameter for calcium indicators. This one may not be used if the
        % selected deconvolution algorithm doesn't need it
        tmp_options.kernel = default_kernel;
        
        % file names for saving avi file
        if save_avi
            tmp_save_avi = sprintf('%sinitialization_%d_%d_%d_%d.avi', log_folder, tmp_block(1), tmp_block(2), tmp_block(3), tmp_block(4));
        else
            tmp_save_avi = save_avi;
        end
        
        % load the patch data
        Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        Ypatch = double(reshape(Ypatch, [], T));
        if nk>1
            Ypatch_dt = detrend_data(Ypatch, nk); % detrend data
            [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch_dt, K, tmp_options, [], tmp_save_avi);
        else
            [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch, K, tmp_options, [], tmp_save_avi);
        end
        
        % put everthing into one struct variable
        tmp_results.center = tmp_center;
        tmp_results.Cn = tmp_Cn;
        tmp_results.PNR = tmp_PNR;
        results{mpatch} = tmp_results;
        %     eval(sprintf('results_patch_%d=tmp_results;', mpatch));  %#ok<PFBFN>
        % display initialization progress
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch);
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
else
    for mpatch=1:(nr_patch*nc_patch) %#ok<*UNRCH>
        % get the indices corresponding to the selected patch
        tmp_patch = patch_pos{mpatch};
        tmp_block = block_pos{mpatch};
        
        % boundaries pixels to be avoided for detecting seed pixels
        tmp_options = options;
        tmp_options.visible_off = true;
        tmp_options.bd = max([(tmp_patch-tmp_block).*[1, -1, 1, -1]; bd, bd, bd, bd], [], 1);
        
        % patch dimension
        tmp_options.d1 = diff(tmp_block(1:2))+1;
        tmp_options.d2 = diff(tmp_block(3:4))+1;
        
        % parameter for calcium indicators. This one may not be used if the
        % selected deconvolution algorithm doesn't need it
        tmp_options.kernel = default_kernel;
        
        % file names for saving avi file
        if save_avi
            tmp_save_avi = sprintf('%sinitialization_%d_%d_%d_%d.avi', log_folder, tmp_block(1), tmp_block(2), tmp_block(3), tmp_block(4));
        else
            tmp_save_avi = save_avi;
        end
        
        % load the patch data
        Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
        Ypatch = double(reshape(Ypatch, [], T));
        if nk>1
            Ypatch_dt = detrend_data(Ypatch, nk); % detrend data
            [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch_dt, K, tmp_options, [], tmp_save_avi);
        else
            [tmp_results, tmp_center, tmp_Cn, tmp_PNR, ~] = greedyROI_endoscope(Ypatch, K, tmp_options, [], tmp_save_avi);
        end
        
        % put everthing into one struct variable
        tmp_results.center = tmp_center;
        tmp_results.Cn = tmp_Cn;
        tmp_results.PNR = tmp_PNR;
        results{mpatch} = tmp_results;
        %     eval(sprintf('results_patch_%d=tmp_results;', mpatch));  %#ok<PFBFN>
        % display initialization progress
        [r, c] = ind2sub([nr_patch, nc_patch], mpatch); 
        fprintf('Patch (%2d, %2d) is done. %2d X %2d patches in total. \n', r, c, nr_patch, nc_patch);
    end
end

%% collect results
for mpatch=1:(nr_patch*nc_patch)
    % get the indices corresponding to the selected patch
    [mr, mc] = ind2sub([nr_patch, nc_patch], mpatch);
    tmp_patch = patch_pos{mr, mc};
    tmp_block = block_pos{mr, mc};
    r0 = tmp_block(1);
    r1 = tmp_block(2);
    c0 = tmp_block(3);
    c1 = tmp_block(4);
    ind_patch = true(r1-r0+1, c1-c0+1);
    ind_patch((tmp_patch(1):tmp_patch(2))-r0+1, (tmp_patch(3):tmp_patch(4))-c0+1) = false;
    
    % unpack results
    tmp_results = results{mpatch};
    %     eval(sprintf('tmp_results=results_patch_%d;', mpatch));
    tmp_Ain = tmp_results.Ain;
    tmp_Ain(ind_patch, :) = 0;
    tmp_Cin = tmp_results.Cin;
    tmp_Cin_raw = tmp_results.Cin_raw;
    tmp_center = tmp_results.center ;
    tmp_Cn = tmp_results.Cn;
    tmp_PNR = tmp_results.PNR;
    if options.deconv_flag
        tmp_Sin = tmp_results.Sin;
        tmp_kernel_pars = tmp_results.kernel_pars;
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
    
    Cn(r0:r1, c0:c1) = max(Cn(r0:r1, c0:c1), tmp_Cn.*(1-ind_patch));
    PNR(r0:r1, c0:c1) = max(PNR(r0:r1, c0:c1), tmp_PNR.*(1-ind_patch));
end
%% export the results
Ain = cell2mat(reshape(Ain, 1, []));
Cin = cell2mat(reshape(Cin, [], 1));
Cin_raw = cell2mat(reshape(Cin_raw, [], 1));
center = cell2mat(reshape(center, [], 1));
obj.A = sparse(Ain);
obj.C = Cin;
obj.C_raw = Cin_raw;
if options.deconv_flag
    Sin = cell2mat(reshape(Sin, [], 1));
    kernel_pars = cell2mat(reshape(kernel_pars, [], 1));
    obj.S = sparse(Sin);
    obj.P.kernel_pars = kernel_pars;
else
    obj.S = zeros(size(obj.C));
end
obj.Cn = Cn;


%% save the results to log
initialization.neuron = obj.obj2struct();
initialization.center = center;
initialization.Cn = Cn;
initialization.PNR = PNR;
log_data.initialization = initialization;

fprintf(flog, '%s\b', get_minute());
fprintf(flog, 'Finished the initialization procedure.\n');
fprintf(flog, 'In total, %d neurons were initialized. \n', size(Ain,2));
fprintf(flog, 'The initialization results were saved as intermediate_results.initialization\n\n');
fclose(flog);