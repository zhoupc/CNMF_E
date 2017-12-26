%% clear the workspace and select data 
clear; clc; close all; 

%% choose multiple datasets or just one  
neuron = Sources2D(); 
nams = {'./data_2p.tif'};          % you can put all file names into a cell array; when it's empty, manually select files 
nams = neuron.select_multiple_files(nams);  %if nam is [], then select data interactively 

%% parameters  
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 8, ...   % GB, memory space you allow to use in MATLAB 
    'memory_size_per_patch', 0.3, ...   % GB, space for loading data within one patch 
    'patch_dims', [30, 40],...  %GB, patch size 
    'batch_frames', 1000);           % number of frames per batch 
   

% -------------------------      SPATIAL      -------------------------  %
gSig = 0.5;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = 10;          % pixel, neuron diameter 
ssub = 1;           % spatial downsampling factor
with_dendrites = false;   % with dendrites or not 
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 20;
    updateA_dist = neuron.options.dist; 
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; %#ok<UNRCH>
    updateA_dist = 5;
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'

% -------------------------      TEMPORAL     -------------------------  %
Fs = 5;             % frame rate
tsub = 1;           % temporal downsampling factor
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -3, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level 
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true);% optimize the baseline); 
nk = 1;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
                    % when changed, try some integers smaller than total_frame/(Fs*30) 
detrend_method = 'local_min';  % compute the local minimum as an estimation of trend. 

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'svd';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
bg_neuron_factor = 0.5;  
ring_radius = round(bg_neuron_factor * gSiz);  % when the ring model used, it is the radius of the ring used in the background model. 
                    %otherwise, it's just the width of the overlapping area 

% -------------------------      MERGING      -------------------------  %
show_merge = false;  % if true, manually verify the merging step 
merge_thr = 0.7;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
dmin = 10;       % minimum distances between two neurons. it is used together with merge_thr
dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only. 

% -------------------------  INITIALIZATION   -------------------------  %
K = [];             % maximum number of neurons per patch. when K=[], take as many as possible 
min_corr = 0.7;     % minimum local correlation for a seeding pixel
min_pnr =10;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 2^2;      % minimum number of nonzero pixels for each neuron
bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames 
save_initialization = false;    % save the initialization procedure as a video. 
use_parallel = true;    % use parallel computation for parallel computing 
show_init = true;   % show initialization results 
choose_params = false; % manually choose parameters 
center_psf = false;  % set the value as true when the background fluctuation is large (usually 1p data) 
                    % set the value as false when the background fluctuation is small (2p)

% ----------------------   MANUAL INTERVENTION  --------------------  %
with_manual_intervention = false; 

% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not 
kt = 1;                 % frame intervals 

% -------------------------    UPDATE ALL    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial -------- 
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'tsub', tsub, ...                       % -------- temporal -------- 
    'deconv_options', deconv_options, ...    '
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background -------- 
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization ----- 
    'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, ...
    'bd', bd, ...
    'center_psf', center_psf);
neuron.Fs = Fs; 

%% distribute data and be ready to run source extraction 
neuron.getReady_batch(pars_envs); 

%% initialize neurons in batch mode 
neuron.initComponents_batch(K, save_initialization, use_parallel); 

%% udpate spatial components for all batches
neuron.update_spatial_batch(use_parallel); 

%% udpate temporal components for all bataches
neuron.update_temporal_batch(use_parallel); 

%% update background 
neuron.update_background_batch(use_parallel); 

%% delete neurons 

%% merge neurons 

%% save workspace 
neuron.save_workspace_batch(); 

%% get the correlation image and PNR image for all neurons 
neuron.correlation_pnr_batch(); 

%% concatenate temporal components 
neuron.concatenate_temporal_batch(); 
neuron.viewNeurons([], neuron.C_raw); 
