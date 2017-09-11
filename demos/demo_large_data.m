clear; clc; close all; 

%% parameters 
% computing environments 
global memory_size_to_use memory_size_per_patch patch_dims ring_radius
memory_size_to_use = 0.1;    % unit: GB. 
memory_size_per_patch = 0.05;  % unit: GB. 
patch_dims = [45, 45];   % unit: GB 
use_gpu = false;    % use GPU for faster matrix factorization 
use_parallel = true;  % use parallel computing 

% neuron, spatial 
gSig = 3; 
gSiz = 11; 
bg_neuron_factor = 1.5; 
ring_radius = round(bg_neuron_factor * gSiz); 
ssub = 1;           % spatial downsampling factor
with_dendrites = false;   % with dendrites or not 

% neuron temporal 
Fs = 6;             % frame rate
tsub = 1;           % temporal downsampling factor
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ... % optimize the baseline
    'optimize_smin', false);  % optimize the threshold 

% merge neurons 
merge_thr = [1e-1, 0.85, 0];     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
center_method = 'mean'; 
dmin = 1; 

%% select data and distribute it into multiple patches. 
% nam = './data_endoscope.tif'; 
nam = './data_endoscope.tif'; 
nam = get_fullname(nam); 

[mat_data, dims] = distribute_data(nam, patch_dims, ring_radius, memory_size_per_patch, memory_size_to_use); 
d1 = dims(1); 
d2 = dims(2); 
numFrames = dims(3); 

neuron = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'gSig', gSig,...    % sigma of the 2D gaussian that approximates cell bodies
    'gSiz', gSiz);      % average neuron size (diameter)
neuron.Fs = Fs;         % frame rate
neuron.file = nam; 
neuron.mat_data = mat_data; 
neuron.options.deconv_options = deconv_options; 

if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron_full.options.search_method = 'dilate'; 
    neuron_full.options.bSiz = 20;
else
    % determine the search locations by selecting a round area
    neuron_full.options.search_method = 'ellipse';
    neuron_full.options.dist = 5;
end

%% initialization 
max_neuron_per_patch = [];  % when [], take as many as possible 
frame_range = [];   % when [], uses all frames 
save_avi = false;    % save the initialization procedure as a video. 

min_corr = 0.8;     % minimum local correlation for a seeding pixel
min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 3^2;      % minimum number of nonzero pixels for each neuron
bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);
neuron.options.nk = 1;  % number of knots for detrending 

neuron.initComponents_parallel(max_neuron_per_patch, frame_range, save_avi); 

