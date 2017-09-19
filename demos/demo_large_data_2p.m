%% clear the workspace and select data 
clear; clc; close all; 

%% choose data 
neuron = Sources2D(); 
nam = []; 
% nam = './data_endoscope.tif'; 
% nam = './msCam1.avi'; 
nam = neuron.select_data(nam);  

%% parameters  
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 8, ...   % GB, memory space you allow to use in MATLAB 
    'memory_size_per_patch', 0.3, ...   % GB, space for loading data within one patch 
    'patch_dims', [64, 64]);  %GB, patch size 
   
% -------------------------      SPATIAL      -------------------------  %
gSig = 1;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = 15;          % pixel, neuron diameter 
ssub = 1;           % spatial downsampling factor
with_dendrites = true;   % with dendrites or not 
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate'; 
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
Fs = 3.7;             % frame rate
tsub = 1;           % temporal downsampling factor
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -3, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level 
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true);% optimize the baseline); 
nk = 10;             % detrending the slow fluctuation. usually 1 is find (no detrending)
                    % when changed, try some integers smaller than total_frame/(Fs*30) 
detrend_method = 'local_min';  % compute the local minimum as an estimation of trend. 

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'svd';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
bg_neuron_factor = 0.5;  
ring_radius = round(bg_neuron_factor * gSiz);  % when the ring model used, it is the radius of the ring used in the background model. 
                    %otherwise, it's just the width of the overlapping area 

% -------------------------      MERGING      -------------------------  %
merge_thr = [1e-1, 0.85, 0];     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'mean';   % method for computing neuron distances 
dmin = 1;  % minimum distances between two neurons 

% -------------------------  INITIALIZATION   -------------------------  %
K = [];             % maximum number of neurons per patch. when K=[], take as many as possible 
min_corr = 0.8;     % minimum local correlation for a seeding pixel
min_pnr = 8;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 1;      % minimum number of nonzero pixels for each neuron
bd = 5;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames 
save_initialization = false;    % save the initialization procedure as a video. 
use_parallel = true;    % use parallel computation for parallel computing 
show_init = true;   % show initialization results 
choose_params = true; % manually choose parameters 
center_psf = false;  % set the value as true when the background fluctuation is large (usually 1p data) 
                    % set the value as false when the background fluctuation is small (2p)
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
neuron.getReady(pars_envs); 

%% initialize neurons from the video data within a selected temporal range 
if choose_params
    % change parameters for optimized initialization
    [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
end

[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel); 
if show_init
    figure;
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
end

%% estimation of background components 
neuron.update_background_parallel(use_parallel); 
% Ybg = neuron.reconstruct_background(); 

%% estimate the temporal components 
neuron.update_temporal_parallel(use_parallel); 

%% estimate the spatial components 
neuron.update_spatial_parallel(use_parallel); 

%% deal with empty neurons 

%% merge neurons 

%% add an ID to each neuron




%% spatial and temporally downsampling 





























