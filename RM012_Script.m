%% Total run script for CNMF_E, for RM012 data
% R. Conor Heins, 05032017

%% clear workspace
clear; clc; close all;  
global numFrames ssub tsub sframe num2read Fs neuron neuron_ds ...
    neuron_full Ybg_weights; % global variables, don't change them manually

%% Configure paths

fprintf('Choose Base CNMF_E Directory...\n')
CNMF_path = uigetdir();
cd(CNMF_path); clear CNMF_path;

run_setup;

%% Choose and load data

% Get user-input to determine whether you're:
% a) loading/processing an existing matfile, which has already been memory-mapped; 
% b) memory-mapping more raw data to an existing matfile, which has
%    already been memory-mapped; 
% c) memory-mapping raw data to a totally new matfile

already_loaded = input('Would you like to load an existing .mat file? (y/n, default n)\n','s');
if isempty(already_loaded)
    already_loaded = 'n';
end

if strcmp(already_loaded,'y')
    fprintf('Choose an existing .mat file to analyze or add data to...\n');
    [fnam, fdir] = uigetfile('*.mat');
    dataFile = fullfile(fdir,fnam);
    data = matfile(dataFile,'Writable',true);
    add_data = input(sprintf('Would you like to add more data to existing file %s ? (y/n, default n)\n',fnam),'s');
    if strcmp(add_data,'y')
        [data,Ysiz,trialsframes] = sequence2mat_mod(data,1);
        origYsiz = data.Ysiz;
        Ysiz_new = [origYsiz(1:2);origYsiz(3)+Ysiz(3)];
        data.Ysiz = Ysiz_new; clear origYsiz Ysiz_new;
        origTrialsFrames = data.trialsframes;
        trialsframes_new = [origTrialsFrames;trialsframes];
        data.trialsframes = trialsframes_new; clear origTrialsFrames trialsframes_new;
    end
else
    fprintf('Choose a folder to determine name of data file...\n');
    dataFile = [uigetdir(),'.mat'];
    [fdir, fnam, temp] = fileparts(dataFile);
    fnam = strcat(fnam,temp); clear temp;
    Y = []; save(dataFile,'Y','-v7.3');
    data = matfile(dataFile,'Writable',true);
    [data,Ysiz,trialsframes] = sequence2mat_mod(data,0);
    data.Ysiz = Ysiz; data.trialsframes = trialsframes;
    add_data = input(sprintf('Would you like to add more data to existing file %s ? (y/n, default n)\n',fnam),'s');
    if strcmp(add_data,'y')
        [data,Ysiz,trialsframes] = sequence2mat_mod(data,1);
        origYsiz = data.Ysiz;
        Ysiz_new = [origYsiz(1:2);origYsiz(3)+Ysiz(3)];
        data.Ysiz = Ysiz_new; clear origYsiz Ysiz_new;
        origTrialsFrames = data.trialsframes;
        trialsframes_new = [origTrialsFrames;trialsframes];
        data.trialsframes = trialsframes_new; clear origTrialsFrames trialsframes_new;
    end
end

%get information about the data
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrames = Ysiz(3);    %total number of frames
trialsframes = data.trialsframes; % trials and frames array 

fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames, and is comprised of %d total trials. \nLoading all data requires %.2f GB RAM\n\n', ...
    d1, d2, numFrames,length(unique(trialsframes(:,1))),prod(Ysiz)*8/(2^30));


%% Create struct with all options and create Source2D class object for storing these and the CNMF-results
ssub = 1; tsub = 1; Fs = 10; 
neuron_full = Sources2D();
neuron_full.Fs = Fs; %frame rate

%Image dimensions and subsampling
options.d1 = d1;
options.d2 = d2;
options.ssub = ssub; %spatial downsampling factor
options.tsub = tsub; %temporal downsampling factor
options.init_method = 'greedy'; % initialization method, either done greedily or with a sparse NMF ('greedy' (default) vs. 'sparse_NMF')

% greedy_corr parameters (greedyROI_corr.m)
options.min_corr = 0.5; %minimum local correlation for initializing a neuron (default: 0.3)
 
% greedyROI parameters (greedyROI.m)
options.gSig        = 6;             % half size of neurons to be found (default: [5,5]) -- also used to compute correlation image
options.gSiz        = 16;            % half size of bounding box for each neuron (default: 2*gSig+1) -- also used to compute correlation image
options.nb          = 1;           % number of background components (default: 1)
options.nIter       = 5;           % maximum number of rank-1 NMF iterations during refining (default: 5)
options.med_app     = 1;           % number of timesteps to be interleaved for fast (approximate) median calculation (default: 1)
options.save_memory = 1;           % process data sequentially to save memory (default: 0)
options.chunkSiz    = 100;         % filter this number of timesteps each time (default: 100)
options.windowSiz   = [32,32];     % size of window over which is computed sequentially (default: 32 x 32)

% sparse_NMF parameters (sparse_NMF_initialization.m)
options.snmf_max_iter = 100;     % max # of sparse NMF iterations (default: 100)
options.err_thr       = 1e-4;    % relative change threshold for stopping sparse_NMF (default: 1e-4)
options.eta           = 1;       % frobenious norm factor: eta*max(Y(:))^2 (default: 1)
options.beta          = 0.8;     % sparsity factor (default: 0.5)

% % HALS parameters (HALS_2d.m)
options.bSiz    = 3;             % expand kernel for HALS growing (default: 3)
options.maxIter = 5;             % maximum number of HALS iterations (default: 5)

% Noise and AR coefficients calculation (preprocess_data.m)
options.noise_range    = [0.25 0.5]; % frequency range over which to estimate the noise (default: [0.25,0.5])
options.noise_method   = 'logmexp';  % method for which to estimate the noise level (default: 'logmexp')
options.flag_g         = false;      % compute global AR coefficients (default: false)
options.lags           = 3;          % number of extra lags when computing the AR coefficients (default: 5)
options.include_noise  = true;          % include early lags when computing AR coefs (default: 0)
options.pixels         = [];         % pixels to include when computing the AR coefs (default: 1:numel(Y)/size(Y,ndims(Y)))
options.split_data     = 1;          % split data into patches for memory reasons (default: 0)
options.block_size     = [64,64];    % block size for estimating noise std in patches (default: [64,64])
options.cluster_pixels = true;       % cluster pixels to active/inactive based on the PSD density (default: true)

% UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
options.search_method  = 'ellipse';      % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'ellipse')
options.use_parallel   = 1;           % update pixels in parallel (default: 1 if present)
 
% determine_search_location.m
options.min_size = 4;                           % minimum size of ellipse axis (default: 3)
options.max_size = 8;                           % maximum size of ellipse axis (default: 8)
options.dist     = 4;                           % expansion factor of ellipse (default: 3)
options.se       = strel('disk',4,0);           % morphological element for dilation (default: strel('disk',4,0))
 
% threshold_components.m
options.nrgthr  = 0.99;                  % energy threshold (default: 0.99)
options.clos_op = strel('square',3);     % morphological element for closing (default: strel('square',3))
options.medw    = [3,3];                 % size of median filter (default: [3,3])

% UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
options.deconv_method = 'constrained_foopsi';     % method for spike deconvolution (default: 'constrained_foopsi')
options.restimate_g   = 1;                        % flag for updating the time constants for each component (default: 1)
options.temporal_iter = 2;                        % number of block-coordinate descent iterations (default: 2)
options.temporal_parallel = 1;                 % flag for parallel updating of temporal components (default: true if present)

% CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
options.method       = 'cvx';   % methods for performing spike inference ('dual','cvx','spgl1','lars') (default:'cvx')
options.bas_nonneg   =     1;   % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 0)
options.fudge_factor =  0.99;   % scaling constant to reduce bias in the time constant estimation (default 0.99 - no scaling)
options.resparse     =     1;   % number of times that the solution is resparsened (default: 0)

% MERGING (merge_ROIs.m)
options.merge_thr  = 0.6;      % merging threshold (default: 0.85)
options.fast_merge = 1;        % flag for using fast merging (default 1)

%  DF/F (extract_DF_F.m)
options.df_prctile = 25;          % percentile to be defined as baseline (default 50, median)
options.df_window  = 100;         % length of running window (default [], no window)

% CONTOUR PLOTS (plot_contours.m)
options.cont_threshold  = 0.9;             %the value above which a specified fraction of energy is explained (default 90%)  

% VIDEO (make_patch_video.m)
options.ind             = [1 2 3 4];       % indeces of components to be shown (deafult: 1:4)
options.skip_frame      = 1;          % skip frames when showing the video (default: 1 (no skipping))
options.sx              = 16;    % half size of representative patches (default: 16)
options.make_avi        = 0;    % flag for saving avi video (default: 0)
options.show_background = 1;    % flag for displaying the background in the denoised panel (default: 1)
options.show_contours   = 1;    % flag for showing the contour plots of the patches in the FoV (default: 0)
options.cmap            = 'default';    % colormap for plotting (default: 'default')
options.name            =['video_',datestr(now,30),'.avi'];    % name of saved video file (default: based on current date)

% PLOT COMPONENTS (view_patches.m)
options.plot_df         = 1;    % flag for displaying DF/F estimates (default: 1)
options.make_gif        = 0;    % save animation (default: 0)
options.save_avi        = 0;    % save video (default: 0)
options.pause_time      = Inf;  % time to pause between each component (default: Inf, user has to click)

% CLASSIFY COMPONENTS (classify components.m)
options.cl_thr          = 0.8;  % overlap threshold for energy for a component to be classified as true (default: 0.8)

% ORDER COMPONENTS (order_components.m)
options.nsd             = 3;    % number of standard deviations (default: 3)
options.nfr             = 5;    % number of consecutive frames (default: 5)

% parameters for microendoscope
options.min_pnr         =  10;    
options.seed_method     = 'auto';  
options.min_pixel       =  5;    % minimum number of nonzero pixels for a neuron
options.bd              =  3;    % number of pixels to be ignored in the boundary
options.deconv_flag     =  1;    % perform deconvolution or not

%Options for running deconvolution
options.deconv_options  = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', false, ... % optimize the baseline
    'optimize_smin', true);  % optimize the threshold 
options.smin            = 5;     % mimimum spike size

neuron_full.options = options;

%% Create convolution kernel to model the shape of calcium transients, and update neuron object accordingly
tau_decay = 1.5;  
tau_rise = 0.1;
nframe_decay = ceil(10*tau_decay*neuron_full.Fs);  % number of frames in decaying period
bound_pars = false;     % bound tau_decay/tau_rise or not
neuron_full.kernel = create_kernel('exp2', [tau_decay, tau_rise]*neuron_full.Fs, nframe_decay, [], [], bound_pars); % add kernel to Sources2D parameters

%% Load the data
sframe=1;						% user input: first frame to read (optional, default:1)
num2read= numFrames;             % user input: how many frames to read (optional, default: until the end)

tic;
if and(ssub==1, tsub==1)
    neuron = neuron_full;
    Y = double(data.Y(:, :, sframe+(1:num2read)-1));
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else
    [Y, neuron_ds] = neuron_full.load_data(nam_mat, sframe, num2read);
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
    neuron = neuron_ds.copy(); 
end
fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);

[d1s, d2s, T] = size(Y); %update d1, d2, and T with new number of frames
Ysiz = [d1s d2s T]';     %update Ysiz with image dimensions

Y = neuron.reshape(Y, 1);       % convert a 3D video into a 2D matrix


%% compute correlation image and peak-to-noise ratio image.
[Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, T, min(T, 500)))));
% show correlation image 
figure('position', [10, 500, 1776, 400]);
subplot(131);
imagesc(Cn, [0, 1]); colorbar;
axis equal off tight;
title('correlation image');

% show peak-to-noise ratio 
subplot(132);
imagesc(pnr,[0,max(pnr(:))*0.98]); colorbar;
axis equal off tight;
title('peak-to-noise ratio');

% show pointwise product of correlation image and peak-to-noise ratio 
subplot(133);
imagesc(Cn.*pnr, [0,max(pnr(:))*0.98]); colorbar;
axis equal off tight;
title('Cn*PNR');

%% Initialization of A,C
debug_on = false;
save_avi = false;
patch_par = [1,1]*1;  % divide the optical field into m X n patches and do initialization patch by patch
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically

min_corr = 0.5;     % minimum local correlation for a seeding pixel
min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 4;      % minimum number of nonzero pixels for each neuron, in terms of the downsampled data
bd = 3;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd, 'deconv_flag', true);
neuron.options.nk = 1;  % number of knots for detrending 

%update 'original' options structure and neuron_raw to reflect new neuron obj's options
options.min_corr = 0.5; options.min_pnr = 10; options.min_pixel = 4; 
options.bd = 3; options.nk = 1; neuron_full.options = options;

% greedy method for initialization
tic;
neuron.options.seed_method = 'auto'; 
[center, Cn, ~] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);

% show results
figure;
imagesc(Cn);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;

% sort neurons
[~, srt] = sort(max(neuron.C, [], 2), 'descend');
neuron.orderROIs(srt);
neuron_init = neuron.copy();

%view and manually delete/trim certain components
% neuron_init.viewNeurons([],neuron_init.C_raw)

%% iteratively update A, C and B
% parameters, merge neurons
display_merge = false;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the background
spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = 1.5;  % spatial range / diameter of neurons

% parameters, estimate the spatial components
update_spatial_method = 'hals_thresh';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
Nspatial = 5;       % this variable has different meanings: 
                    %1) udpate_spatial_method=='hals' or 'hals_thresh',
                    %then Nspatial is the maximum iteration 
                    %2) update_spatial_method== 'nnls', it is the maximum
                    %number of neurons overlapping at one pixel 
                    
neuron.options.maxIter = 2;   % iterations to update C

% parameters for running iteratiosn 
nC = size(neuron.C, 1);    % number of neurons 

maxIter = 3;        % maximum number of iterations 
miter = 1; 
while miter <= maxIter
    %% merge neurons, order neurons and delete some low quality neurons
     if miter ==1
        merge_thr = [1e-1, 0.8, .1];     % thresholds for merging neurons
        % corresponding to {sptial overlaps, temporal correlation of C,
        %temporal correlation of S}
    else
        merge_thr = [0.6, 0.5, 0.1]; 
     end
    
    % merge neurons
    cnmfe_quick_merge;              % run neuron merges
    
    %% udpate background (cell 1, the following three blocks can be run iteratively)
    % estimate the background
    tic;
    
    %ALSO KNOWN AS: 'cnmfe_update_BG'
    Ybg = Y-neuron.A*neuron.C;
    rr = ceil(neuron.options.gSiz * bg_neuron_ratio);
    active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
    [Ybg, Ybg_weights] = neuron.localBG(Ybg, spatial_ds_factor, rr, active_px, neuron.P.sn, thresh); % estiamte local background.
    % subtract the background from the raw data.
    Ysignal = Y - Ybg;
    
    % estimate noise
    if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
        %% estimate the noise for all pixels
        b0 =zeros(size(Ysignal,1), 1);
        sn = b0;
        parfor m=1:size(neuron.A,1)
            [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m, :));
        end
        Ysigma = bsxfun(@minus, Ysignal, b0);
        neuron.P.sn = sn;
    end
    
    fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
    % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.
    
    %% update spatial & temporal components
    tic;
    for m=1:5    
        %temporal
        neuron.updateTemporal_endoscope(Ysignal);
        cnmfe_quick_merge;              % run neuron merges
        %spatial
        neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
        neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        if isempty(merged_ROI)
            break;
        end
    end
    fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
    
    %% pick neurons from the residual (cell 4).
    if miter==1
        
        %reset thresholds for picking neurons from residual
        min_corr = 0.6;     % minimum local correlation for a seeding pixel
        min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
        min_pixel = 6;      % minimum number of nonzero pixels for each neuron
        bd = 10;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
        neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
        'min_pixel', min_pixel, 'bd', bd, 'deconv_flag', true);
    
        neuron.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
        [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par, 'auto'); % method can be either 'auto' or 'manual'
    end
    
    %% stop the iteration 
    temp = size(neuron.C, 1); 
    if or(nC==temp, miter==maxIter)
        break; 
    else
        miter = miter+1; 
        nC = temp; 
    end
end

%% view neurons
if view_neurons
    neuron.viewNeurons([], neuron.C_raw);
end

% RUN ON FULL RESOLUTION AFTER DOWNSAMPLING
% ALSO KNOWN AS: 'cnmfe_full'
%% upsample the CNMF results 

if or(ssub>1, tsub>1)
    neuron_ds = neuron.copy();  % save the result
    neuron = neuron_full.copy();
    A0 = neuron_ds.reshape(neuron_ds.A, 2);
    C0 = neuron_ds.C;
    C0_raw = neuron_ds.C_raw;
    kernel_pars = neuron_ds.P.kernel_pars;
    K = size(C0, 1);     % number of neurons

    neuron.A = neuron.reshape(imresize(A0, [d1, d2]), 1);
    C = zeros(K, num2read);
    C(:, 1:T*tsub) = resample(C0', tsub, 1)';
    temp = num2read - T*tsub;
    if temp>0
        C(:, (num2read-temp+1):end) = C(:, T*tsub) * ones(1, temp);
    end
    neuron.C = C;
    neuron.C_raw = zeros(K, num2read);
    neuron.S = zeros(K, num2read);
    neuron.P.kernel_pars = kernel_pars * tsub;
    
    %% load data
    Y = data.Y(:, :, sframe:(sframe+num2read-1));
    Y = neuron.reshape(double(Y), 1);
    
    if ssub ==1
        neuron.P.sn = sn;
    elseif ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
        sn = neuron.estNoise(Y);
        neuron.P.sn = sn;
    else
        sn = neuron.P.sn;
    end
    
    %% estimate the background
    %ALSO KNOWN AS: 'cnmfe_update_BG'
    tic
    
    Ybg = Y-neuron.A*neuron.C;
    rr = ceil(neuron.options.gSiz * bg_neuron_ratio);
    active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
    [Ybg, Ybg_weights] = neuron.localBG(Ybg, spatial_ds_factor, rr, active_px, neuron.P.sn, thresh); % estiamte local background.
    % subtract the background from the raw data.
    Ysignal = Y - Ybg;
    
    % estimate noise
    if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
        %% estimate the noise for all pixels
        b0 =zeros(size(Ysignal,1), 1);
        sn = b0;
        parfor m=1:size(neuron.A,1)
            [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m, :));
        end
        Ysigma = bsxfun(@minus, Ysignal, b0);
        neuron.P.sn = sn;
    end
    
     fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
     
    %% update the spatial and components
    maxIter = 1;
    for miter=1:maxIter
        fprintf('Iteration %d/%d to update spatial and temporal components\n', miter, maxIter);
        tic;
        neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
        neuron.trimSpatial(0.02);
        fprintf('Time cost in updating neuronal spatial components:     %.2f seconds\n', toc);
        
        neuron.updateTemporal_endoscope(Ysignal);
        fprintf('Time cost in updating neuronal temporal components:     %.2f seconds\n', toc);
    end

neuron_full = neuron.copy();
end

%% display contours of the neurons
figure;
Cnn = correlation_image(neuron_full.reshape(Ysignal(:, 1:5:end), 2), 4);
neuron.Coor = plot_contours(neuron.A, Cnn, 0.8, 0, [], [], 2);
colormap winter;
axis equal; axis off;
title('contours of estimated neurons');

%% view and edit (if you want) neurons one by one
if viewNeurons
    neuron_full.viewNeurons([],neuron_full.C_raw);
end
%% plot centers of neurons

centers = neuron_full.estCenter();

figure;
imagesc(Cnn);
hold on; scatter(centers(:,2),centers(:,1),7,'or','filled');
colormap; axis off tight equal;

%% quickly merge neighbors that are too close to eachother
[merged_ROIs, newIDs] = neuron_full.MergeNeighbors(5, 'maximum');
new_centers = neuron_full.estCenter();

figure;
imagesc(Cnn);
hold on; scatter(new_centers(:,2),new_centers(:,1),7,'or','filled');
colormap; axis off tight equal;


%% save results

saveName = [dir_nm '_results.mat'];
save(saveName,'neuron_full','d1','d2','numFrames','ssub','tsub','trials_frames');

