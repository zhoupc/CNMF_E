%% clear the workspace and select data 
clear; clc; close all; 

%% choose data 
neuron = Sources2D(); 
% nam = []; 
% nam = './data_endoscope.tif'; 
nam = './msCam1.avi'; 
nam = neuron.select_data(nam);  

%% parameters  
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 4, ...   % GB, memory space you allow to use in MATLAB 
    'memory_size_per_patch', 0.5, ...   % GB, space for loading data within one patch 
    'patch_dims', [45, 45]);  %GB, patch size 
   
% -------------------------      SPATIAL      -------------------------  %
gSig = 4;           % pixel, gaussian width of a gaussian kernel that approximating a typical neuron 
gSiz = 17;          % pixel, neuron diameter 
ssub = 1;           % spatial downsampling factor
with_dendrites = false;   % with dendrites or not 
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 20;
    updateA_dist = neuron.options.dist; 
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse';
    updateA_dist = 5;
    updateA_bSiz = neuron.options.dist;
end

% -------------------------      TEMPORAL     -------------------------  %
Fs = 6;             % frame rate
tsub = 1;           % temporal downsampling factor
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ... % optimize the baseline
    'optimize_smin', false);  % optimize the threshold 

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
bg_neuron_factor = 1.5;  
ring_radius = round(bg_neuron_factor * gSiz);  % radius of the ring used in the background model 

% -------------------------      MERGING      -------------------------  %
merge_thr = [1e-1, 0.85, 0];     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'mean';   % method for computing neuron distances 
dmin = 1;  % minimum distances between two neurons 

% -------------------------  INITIALIZATION   -------------------------  %
K = [];             % maximum number of neurons per patch. when K=[], take as many as possible 
min_corr = 0.7;     % minimum local correlation for a seeding pixel
min_pnr = 5;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 3^2;      % minimum number of nonzero pixels for each neuron
bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames 
save_initialization = false;    % save the initialization procedure as a video. 
show_init = true;   % show initialization results 
% -------------------------    UPDATE ALL    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial -------- 
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'tsub', tsub, ...                       % -------- temporal -------- 
    'deconv_options', deconv_options, ...    
    'background_model', bg_model, ...       % -------- background -------- 
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization ----- 
    'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, ...
    'bd', bd);

%% distribute data and be ready to run source extraction 
neuron.getReady(pars_envs); 


%% initialize neurons from the video data within a selected temporal range 
[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization); 
if show_init
    figure;
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 5);
end

%% estimation of background components 



