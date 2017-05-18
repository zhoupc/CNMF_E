% Total run script for CNMF_E, for RM012 data
% R. Conor Heins, 05032017

%% clear workspace
clear; clc; close all;  
global numFrames ssub tsub sframe num2read Fs ...
    neuron_full; % global variables, don't change them manually

%% Configure paths

fprintf('Choose Base CNMF_E Directory...\n')
CNMF_path = uigetdir();
cd(CNMF_path); clear CNMF_path;

cnmfe_setup;

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

process_patches = 1; %flag for whether you want to run CNMF_E on individual spatial patches of data (can be done in parallel)
patch_size = [150, 150]; %patch size
overlap = [25 25]; %patch overlap
min_patch_sz = [26,26];
patches = construct_patches(Ysiz(1:end-1),patch_size,overlap,min_patch_sz);

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
options.min_corr = 0.3; %minimum local correlation for initializing a neuron (default: 0.3)
options.nk       = 1;   %number of knots in B-Spline basis for detrending
 
% greedyROI parameters (greedyROI.m)
options.gSig        = 6;             % half size of neurons to be found (default: [5,5]) -- also used to compute correlation image
options.gSiz        = 16;            % half size of bounding box for each neuron (default: 2*gSig+1) -- also used to compute correlation image
options.nb          = 1;           % number of background components (default: 1)
options.nIter       = 5;           % maximum number of rank-1 NMF iterations during refining (default: 5)
options.med_app     = 1;           % number of timesteps to be interleaved for fast (approximate) median calculation (default: 1)
options.save_memory = 0;           % process data sequentially to save memory (default: 0)
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
options.lags           = 5;          % number of extra lags when computing the AR coefficients (default: 5)
options.include_noise  = true;          % include early lags when computing AR coefs (default: 0)
options.pixels         = [];         % pixels to include when computing the AR coefs (default: 1:numel(Y)/size(Y,ndims(Y)))
options.split_data     = 0;          % split data into patches for memory reasons (default: 0)
options.block_size     = [64,64];    % block size for estimating noise std in patches (default: [64,64])
options.cluster_pixels = true;       % cluster pixels to active/inactive based on the PSD density (default: true)

% UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
options.search_method  = 'ellipse';      % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'ellipse')
options.use_parallel   = 0;           % update pixels in parallel (default: 1 if present)
 
% determine_search_location.m
options.min_size = 3;                           % minimum size of ellipse axis (default: 3)
options.max_size = 6;                           % maximum size of ellipse axis (default: 8)
options.dist     = 2;                           % expansion factor of ellipse (default: 3)
options.se       = strel('disk',3,0);           % morphological element for dilation (default: strel('disk',4,0))
 
% threshold_components.m
options.nrgthr  = 0.99;                  % energy threshold (default: 0.99)
options.clos_op = strel('square',3);     % morphological element for closing (default: strel('square',3))
options.medw    = [3,3];                 % size of median filter (default: [3,3])

% UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
options.deconv_method = 'constrained_foopsi';     % method for spike deconvolution (default: 'constrained_foopsi')
options.restimate_g   = 1;                        % flag for updating the time constants for each component (default: 1)
options.temporal_iter = 2;                        % number of block-coordinate descent iterations (default: 2)
options.temporal_parallel = false;                 % flag for parallel updating of temporal components (default: true if present)

% CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
options.method       = 'cvx';   % methods for performing spike inference ('dual','cvx','spgl1','lars') (default:'cvx')
options.bas_nonneg   =     1;   % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 0)
options.fudge_factor =  0.99;   % scaling constant to reduce bias in the time constant estimation (default 0.99 - no scaling)
options.resparse     =     1;   % number of times that the solution is resparsened (default: 0)

% MERGING (merge_ROIs.m)
options.merge_thr  = 0.85;      % merging threshold (default: 0.85)
options.fast_merge = 1;        % flag for using fast merging (default 1)

%  DF/F (extract_DF_F.m)
options.df_prctile = 25;          % percentile to be defined as baseline (default 50, median)
options.df_window  = 100;         % length of running window (default [], no window)

% CONTOUR PLOTS (plot_contours.m)
options.cont_threshold  = 0.9;             %the value above which a specified fraction of energy is explained (default 90%)  

% VIDEO (make_patch_video.m)
options.ind             = [1 2 3 4];       % indeces of components to be shown (deafult: 1:4)
options.skip_frame      = 1;          % skip frames when showing the video (default: 1 (no skipping))
options.sx              = 18;    % half size of representative patches (default: 16)
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
options.min_pixel       =  8;    % minimum number of nonzero pixels for a neuron
options.bd              =  0;    % number of pixels to be ignored in the boundary
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
tau_decay = 1;  
tau_rise = 0.1;
nframe_decay = ceil(10*tau_decay*neuron_full.Fs);  % number of frames in decaying period
bound_pars = false;     % bound tau_decay/tau_rise or not
neuron_full.kernel = create_kernel('exp2', [tau_decay, tau_rise]*neuron_full.Fs, nframe_decay, [], [], bound_pars); % add kernel to Sources2D parameters

%% Load and analyze the data 
sframe=1;						% user input: first frame to read (optional, default:1)
num2read= numFrames;             % user input: how many frames to read (optional, default: until the end)

if process_patches
    
    RESULTS(length(patches)) = struct();
    
    %% PARALLEL CNMF_E
    
    parfor i = 1:length(patches)
        
        neuron_slave = neuron_full.copy();
        sframe_patch = sframe; num2read_patch = num2read;
        Y = double(data.Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),sframe_patch+(1:num2read_patch)-1));
        [d1,d2,T] = size(Y);
        neuron_slave.options.d1 = d1; neuron_slave.options.d2 = d2;
        Y = neuron_slave.reshape(Y, 1);       % convert a 3D video into a 2D matrix
        
        %% Initialization of A,C -- parameters
        patch_par = [1,1]*1;
        K = [];
        
        %% iteratively update A, C and B -- parameters
        % parameters, merge neurons
        
        % parameters, estimate the background
        spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
        thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
        
        bg_neuron_ratio = 2;  % spatial range / diameter of neurons
        
        % parameters, estimate the spatial components
        update_spatial_method = 'hals_thresh';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
        Nspatial = 5;       % this variable has different meanings:
        %1) udpate_spatial_method=='hals' or 'hals_thresh',
        %then Nspatial is the maximum iteration
        %2) update_spatial_method== 'nnls', it is the maximum
        %number of neurons overlapping at one pixel
        
        neuron_slave.options.maxIter = 2;   % iterations to update C
        
        %% Initialize spatial and temporal components
        
        neuron_slave.initComponents_endoscope(Y, K, patch_par, [], []);
        
        %% Continue with full CNMF if neurons are found in patch
        
        if size(neuron_slave.A,2)>0
            [~, srt] = sort(max(neuron_slave.C, [], 2), 'descend');
            neuron_slave.orderROIs(srt);
            
            neuron_slave.options.maxIter = 3;   % iterations to update C
            
            % parameters for running iteratiosn
            nC = size(neuron_slave.C, 1);    % number of neurons
            
            maxIter = 3;        % maximum number of iterations
            miter = 1;
            while miter <= maxIter
                %% merge neurons, order neurons and delete some low quality neurons
                if miter ==1
                    merge_thr = [1e-1, 0.8, .1];     % thresholds for merging neurons
                    % corresponding to {sptial overlaps, temporal correlation of C,
                    %temporal correlation of S}
                else
                    merge_thr = [0.4, 0.5, 0.1];
                end
                
                if nC <= 1
                    break
                end
                
                % merge neurons
                neuron_slave.quickMerge(merge_thr); % run neuron merges
                %sort neurons
                [~,srt] = sort(max(neuron_slave.C,[],2).*max(neuron_slave.A,[],1)','descend');
                neuron_slave.orderROIs(srt);
                
                %% udpate background (cell 1, the following three blocks can be run iteratively)
                
                Ybg = Y-neuron_slave.A*neuron_slave.C;
                rr = ceil(neuron_slave.options.gSiz * bg_neuron_ratio);
                active_px = []; %(sum(IND, 2)>0); %If some missing neurons are not covered by active_px, use [] to replace IND
                [Ybg, Ybg_weights] = neuron_slave.localBG(Ybg,spatial_ds_factor,rr,active_px,neuron_slave.P.sn,thresh); %estimate local background
                
                %subtract background from the raw data to obtain signal for
                %subsequent CNMF
                Ysignal = Y - Ybg;
            
                % estimate noise
                if ~isfield(neuron_slave.P,'sn') || isempty(neuron_slave.P.sn)
                    % estimate the noise for all pixels
                    b0 = zeros(size(Ysignal,1), 1);
                    sn = b0;
                    for m = 1:size(neuron_slave.A,1)
                        [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m,:));
                    end
                    Ysigma = bsxfun(@minus,Ysignal,b0);
                end
                
                %% update spatial & temporal components
                tic;
                for m=1:2
                    %temporal
                    neuron_slave.updateTemporal_endoscope(Ysignal);
                    
                    % merge neurons
                    [merged_ROI, ~] = neuron_slave.quickMerge(merge_thr); % run neuron merges
                    %sort neurons
                    [~,srt] = sort(max(neuron_slave.C,[],2).*max(neuron_slave.A,[],1)','descend');
                    neuron_slave.orderROIs(srt);
                    
                    %spatial
                    neuron_slave.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
                    neuron_slave.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
                    if isempty(merged_ROI)
                        break;
                    end
                end
                fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
                
                %% pick neurons from the residual (cell 4).
                if miter==1
                    
                    %reset thresholds for picking neurons from residual
                    min_corr = 0.8;     % minimum local correlation for a seeding pixel
                    min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
                    min_pixel = 8;      % minimum number of nonzero pixels for each neuron
                    bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
                    neuron_slave.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
                        'min_pixel', min_pixel, 'bd', bd, 'deconv_flag', true);
                    
                    neuron_slave.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
                    [center_new, Cn_res, pnr_res] = neuron_slave.pickNeurons(Ysignal - neuron_slave.A*neuron_slave.C, patch_par, 'auto'); % method can be either 'auto' or 'manual'
                end
                
                %% stop the iteration
                temp = size(neuron_slave.C, 1);
                if or(nC==temp, miter==maxIter)
                    break;
                else
                    miter = miter+1;
                    nC = temp;
                end
            end
        end
        
        RESULTS(i).A = neuron_slave.A;
        RESULTS(i).C = neuron_slave.C;
        RESULTS(i).C_raw = neuron_slave.C_raw; 
        RESULTS(i).S = neuron_slave.S;
        RESULTS(i).P = neuron_slave.P;
        fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(length(patches)), '.\n']);
        
    end
    
    %% STITCH RESULTS FROM EACH PATCH BACK TOGETHER
    neuron_full.P.kernel_pars = [];
    for i = 1:length(patches)
        if size(RESULTS(i).A,2)>0
            Achunk = zeros(neuron_full.options.d1,neuron_full.options.d2,size(RESULTS(i).A,2));
            for k = 1:size(RESULTS(i).A,2) 
                Achunk(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),k)=reshape(RESULTS(i).A(:,k),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
            end
             neuron_full.A = [neuron_full.A,reshape(Achunk,d1*d2,k)];
             neuron_full.C = [neuron_full.C;RESULTS(i).C];
             neuron_full.C_raw = [neuron_full.C_raw;RESULTS(i).C_raw];
             neuron_full.S = [neuron_full.S;RESULTS(i).S];
             neuron_full.P.sn = [neuron_full.P.sn,RESULTS(i).P.sn];
             neuron_full.P.kernel_pars = [neuron_full.P.kernel_pars;RESULTS(i).P.kernel_pars];      
        end
        clear Achunk;
    end
    
    
end

%% merge neurons with new merge_thr criterion to merge overlapping neurons at patch-intersections

display_merge = 1;

display_merge = 1;
merge_thr = [0.3 0.3 0.1];
neuron_bk = neuron_full.copy();
[merged_ROI, newIDs] = neuron_full.quickMerge(merge_thr);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
% A: spatial shapes; S: spike counts; C: calcium traces
if display_merge && ~isempty(merged_ROI)
    figure('position', [1,1, 1200, 600]);
    ind_before = false(size(neuron_bk.A, 2), 1);
    ind_after = false(size(neuron_full.A, 2), 1);
    m = 1;
    while m<=length(merged_ROI)
        subplot(221);
        tmp_img = neuron_bk.overlapA(merged_ROI{m});
        imagesc(tmp_img);
        axis equal off tight;
        subplot(222);
        imagesc(tmp_img);
        axis equal off tight;
        [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
        xlim([min(tmp_c)-10, max(tmp_c)+10]);
        ylim([min(tmp_r)-10, max(tmp_r)+10]);
        axis off;
        subplot(2,2,3:4);
        tmp_C = neuron_bk.C_raw(merged_ROI{m}, :)';
        plot(tmp_C, 'linewidth', 2);
        
        temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
        if strcmpi(temp, 'n')
            ind_after(newIDs(m)) = true;
            ind_before(merged_ROI{m}) = true;
            m = m+1;
        elseif strcmpi(temp, 'b')
            m = m-1;
        elseif strcmpi(temp, 'e')
            break;
        else
            m = m+1;
        end
    end
    
    neuron_full.A = [neuron_full.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
    neuron_full.C = [neuron_full.C(~ind_after, :); neuron_bk.C(ind_before, :)];
    neuron_full.C_raw = [neuron_full.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
    neuron_full.S = [neuron_full.S(~ind_after, :); neuron_bk.S(ind_before, :)];
    neuron_full.P.kernel_pars = [neuron_full.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
    close; 
end

% sort neurons
[Cpnr, srt] = sort(max(neuron_full.C, [], 2).*max(neuron_full.A, [], 1)', 'descend');
neuron_full.orderROIs(srt);

%% view neurons

view_neurons = 1;

if view_neurons
    neuron_full.viewNeurons([], neuron_full.C_raw);
end

%% compute correlation image and peak-to-noise ratio image with centers of neurons

centers = neuron_full.estCenter();

T = 5000; %end frame for computing correlation image
Y = double(data.Y(:,:,1:T));
Y = neuron_full.reshape(Y,1);

[Cn, pnr] = neuron_full.correlation_pnr(Y(:, round(linspace(1, T, min(T, 500)))));
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
hold on; scatter(centers(:,2),centers(:,1),7,'or','filled');
axis equal off tight;
title('Cn*PNR');

neuron_full.Cn = Cn;


%% save results

saveName = [fdir,filesep,fnam(1:end-4),'_results.mat'];
save(saveName,'neuron_full','patches','numFrames','ssub','tsub','trialsframes');