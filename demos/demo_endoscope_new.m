%% clear workspace
clear; clc; close all;
global  d1 d2 numFrame ssub tsub sframe num2read Fs neuron neuron_ds ...
    neuron_full Ybg_weights; %#ok<NUSED> % global variables, don't change them manually

%% select data and map it to the RAM
% nam = '~/Dropbox/github/CNMF_E/demos/data_endoscope.tif';
cnmfe_choose_data;

%% create Source2D class object for storing results and parameters
Fs = 6;             % frame rate
ssub = 1;           % spatial downsampling factor
tsub = 1;           % temporal downsampling factor
gSig = 4;           %width of the gaussian kernel, which can approximates the average neuron shape
gSiz = 15;          % maximum diameter of neurons in the image plane. larger values are preferred.
neuron_full = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'ssub', ssub, 'tsub', tsub, ...  % downsampleing
    'gSig', gSig,...
    'gSiz', gSiz);
neuron_full.Fs = Fs;         % frame rate

% create convolution kernel to model the shape of calcium transients
tau_decay = 1;  %
tau_rise = 0.1;
nframe_decay = ceil(6*tau_decay*neuron_full.Fs);  % number of frames in decaying period
bound_pars = false;     % bound tau_decay/tau_rise or not
neuron_full.kernel = create_kernel('exp2', [tau_decay, tau_rise]*neuron_full.Fs, nframe_decay, [], [], bound_pars);

%% downsample data for fast and better initialization
sframe=1;						% user input: first frame to read (optional, default:1)
num2read= numFrame;             % user input: how many frames to read   (optional, default: until the end)

tic;
cnmfe_load_data;
fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);

Y = neuron.reshape(Y, 1);       % convert a 3D video into a 2D matrix

%% compute correlation image and peak-to-noise ratio image.
cnmfe_show_corr_pnr;    % this step is not necessary, but it can give you some...
                        % hints on parameter selection, e.g., min_corr & min_pnr

%% initialization of A, C
% parameters
debug_on = false;
save_avi = false;
patch_par = [2,2]; %1;  % divide the optical field into m X n patches and do initialization patch by patch
K = 300; % maximum number of neurons to search within each patch. you can use [] to search the number automatically

min_corr = 0.85;     % minimum local correlation for a seeding pixel
min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 5;      % minimum number of nonzero pixels for each neuron
bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);

% greedy method for initialization
tic;
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);

% show results
figure;
imagesc(Cn);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;

% sort neurons
[~, srt] = sort(max(neuron.C, [], 2)./get_noise_fft(neuron.C), 'descend');
neuron.orderROIs(srt);
neuron_init = neuron.copy();

%% iteratively update A, C and B
% parameters, merge neurons
display_merge = false;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the background
spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
thresh = 5;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
    sn = neuron.estNoise(Y);
else
    sn = neuron.P.sn; 
end
bg_neuron_ratio = 2;  % spatial range / diameter of neurons

% parameters, estimate the spatial components
maxIter_spatial = 10;       % number of iterations required
with_dendrites = false;
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron.options.search_method = 'dilate';  %#ok<UNRCH>
    neuron.options.bSiz = 5;
else
    % determine the search locations by selecting a round area
    neuron.options.search_method = 'ellipse';
    neuron.options.dist = 3;
end

% parameters, estimate the temporal components
smin = 4;       % thresholding the amplitude of the spike counts as smin*noise level

neuron.options.maxIter = 4;   % iterations to update C

% parameters for running iteratiosn 
nC = size(neuron.C, 2);    % number of neurons 

maxIter = 5;        % maximum number of iterations 
miter = 1; 
while miter <= maxIter
    %% merge neurons, order neurons and delete some low quality neurons
    % parameters
    if miter<2
        merge_thr = [0, 0.7, 0];     % thresholds for merging neurons
        % corresponding to {sptial overlaps, temporal correlation of C,
        %temporal correlation of S}
    else
        merge_thr = [.1, 0.6, 0];
    end
    
    % merge neurons
    cnmfe_quick_merge;              % run neuron merges
    
    %% udpate background (cell 1, the following three blocks can be run iteratively)
    % estimate the background
    tic;
    cnmfe_update_BG;
    fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
    % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.
    
    %% update spatial components
    tic;
    neuron.updateSpatial_endoscope(Ysignal, maxIter_spatial);
    fprintf('Time cost in updating neuronal spatial components:     %.2f seconds\n', toc);
    
    %% update temporal components.
    tic;
    neuron.updateTemporal_endoscope(Ysignal, smin);
    fprintf('Time cost in updating neuronal temporal components:     %.2f seconds\n', toc);
    
    %% pick neurons from the residual (cell 4).
    if miter==1
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

%% apply results to the full resolution
if or(ssub>1, tsub>1)
    neuron_ds = neuron.copy();  % save the result
    neuron = neuron_full.copy();
    cnmfe_full;
    neuron_full = neuron.copy();
end

%% display neurons
dir_neurons = sprintf('%s%s%s_neurons%s', dir_nm, filesep, file_nm, filesep);
if exist('dir_neurons', 'dir')
    temp = cd();
    cd(dir_neurons);
    delete *;
    cd(temp);
else
    mkdir(dir_neurons);
end
neuron.viewNeurons([], neuron.C_raw, dir_neurons);
close(gcf); 

%% display contours of the neurons
figure;
Cnn = correlation_image(Ysignal(:, 1:5:end), 4, d1, d2);
neuron.viewContours(Cnn, 0.9, 0); % correlation image computed with background-subtracted data
colormap winter;
axis equal; axis off;
title('contours of estimated neurons');

% plot contours with IDs
[Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, T, min(T, 1000)))));
figure;
plot_contours(neuron.A, Cn, 0.9, 1, [], neuron.Coor, 2);
colormap winter;
title('contours of estimated neurons');

%% check spatial and temporal components by playing movies
% save_avi = false;
% avi_name = 'play_movie.avi';
% neuron.Cn = Cn;
% neuron.runMovie(Ysignal, [0, 50], save_avi, avi_name);

%% save video
kt = 3;     % play one frame in every kt frames
save_avi = false;

cnmfe_save_video;

%% save results
globalVars = who('global');
eval(sprintf('save %s%s%s_results.mat %s', dir_nm, filesep, file_nm, strjoin(globalVars)));