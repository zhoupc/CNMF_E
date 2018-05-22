%% Compile singing/sleep data in one day
Datafolder='/net/feevault/data0/elm/ProcessedCalciumData/sleep/6962/042917/';
singing = matfile(fullfile(Datafolder,'compiled.mat'));
sleep=matfile(fullfile(Datafolder,'CompiledSleep_6962_20170429_105506-20170430_104328.mat'));
Myfolder='/home/shijiegu/SingAndSleep';


Bird=sleep.Bird;
SleepStart=sleep.SleepStart;
Day=datestr(SleepStart(1:13),'ddmmmyyyy');

Y = cat(3,singing.Y,sleep.Y);
display('Y cated')
VIDEOfs=singing.VIDEOfs;
FILE = fullfile(Myfolder, [Bird '_' 'SingANDSleep_' Day]); 
save(FILE, 'VIDEOfs','Y','-v7.3');

%% Demo Starts
% parameter sweep version 
global  d1 d2 numFrame ssub tsub sframe num2read Fs neuron neuron_ds ...
    neuron_full Ybg_weights outputfolder data Ysiz nam min_pnr gSiz gSig; %#ok<NUSED> % global variables, don't change them manually

addpath(genpath('/home/shijiegu/CNMF_E-master'))

%% select data and map it to the RAM
%nam = '/Users/gushijie/Documents/MATLAB/CalciumExtraction/CaELM_row3888.mat';
nam=FILE;
%nam='/Volumes/data0-shared/elm/ProcessedCalciumData/6865_Jan9/compiled.mat';
%cnmfe_choose_data;
[dir_nm, file_nm, file_type] = fileparts(nam);
nam_mat = nam;

% information of the data 
data = matfile(nam_mat);
Ysiz = size(data.Y); 
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames

fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));
%% create Source2D class object for storing results and parameters
Fs = 30;             % frame rate
ssub = 1;           % spatial downsampling factor
tsub = 1;           % temporal downsampling factor
gSig = gSig;        % width of the gaussian kernel, which can approximates the average neuron shape
gSiz = gSiz;        % maximum diameter of neurons in the image plane. larger values are preferred.
neuron_full = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'ssub', ssub, 'tsub', tsub, ...  % downsampleing
    'gSig', gSig,...    % sigma of the 2D gaussian that approximates cell bodies
    'gSiz', gSiz);      % average neuron size (diameter)
neuron_full.Fs = Fs;         % frame rate

% with dendrites or not 
with_dendrites = false;
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron_full.options.search_method = 'dilate'; 
    neuron_full.options.bSiz = 20;
else
    % determine the search locations by selecting a round area
    neuron_full.options.search_method = 'ellipse';
    neuron_full.options.dist = 5;
end

%% options for running deconvolution 
neuron_full.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', false, ... % optimize the baseline
    'optimize_smin', true);  % optimize the threshold 

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
debug_on = false;   % visualize the initialization procedue. 
save_avi = false;   %save the initialization procedure as an avi movie. 
patch_par = [1,1]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large 
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically

min_corr = 0.85;     % minimum local correlation for a seeding pixel
min_pnr = min_pnr;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 5;      % minimum number of nonzero pixels for each neuron
bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);
neuron.options.nk = 1;  % number of knots for detrending 

% greedy method for initialization
tic;
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);

% show results
% figure;
% imagesc(Cn, [0.1, 0.95]);
% hold on; plot(center(:, 2), center(:, 1), 'or');
% colormap; axis off tight equal;

% sort neurons
[~, srt] = sort(max(neuron.C, [], 2), 'descend');
neuron.orderROIs(srt);
neuron_init = neuron.copy();
display('line 99')
%% iteratively update A, C and B
% parameters, merge neurons
display_merge = false;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the background
spatial_ds_factor = 1;      % spatial downsampling factor. it's for faster estimation
thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = 1;  % spatial range / diameter of neurons

% parameters, estimate the spatial components
update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
Nspatial = 5;       % this variable has different meanings: 
                    %1) udpate_spatial_method=='hals' or 'hals_thresh',
                    %then Nspatial is the maximum iteration 
                    %2) update_spatial_method== 'nnls', it is the maximum
                    %number of neurons overlapping at one pixel 
               
% parameters for running iteratiosn 
nC = size(neuron.C, 1);    % number of neurons 

maxIter = 2;        % maximum number of iterations 
miter = 1; 
display('line 114')
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
    cnmfe_update_BG;
    fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
    % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.
    
    %% update spatial & temporal components
    tic;
    for m=1:2  
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
        neuron.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
        [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par); % method can be either 'auto' or 'manual'
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


%% delete some neurons and run CNMF-E iteration 
% neuron.viewNeurons([], neuron.C_raw); 
tic;
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
%update spatial & temporal components
tic;
for m=1:2
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
    cnmfe_quick_merge;              % run neuron merges
    %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% display neurons

teststring=sprintf('%s_%s_singsleepresults.mat', Bird, Day);
neuron.viewNeurons([], neuron.C_raw, teststring);
close(gcf);

ColorAllNeurons(neuron.A)
fprintf('Almost done');

% %% display contours of the neurons
% neuron.Coor = neuron.get_contours(0.8); % energy within the contour is 80% of the total 
% figure;
% Cnn = correlation_image(neuron.reshape(Ysignal(:, 1:5:end), 2), 4);
% neuron.Coor = plot_contours(neuron.A, Cnn, 0.8, 0, [], neuron.Coor, 2);
% colormap winter;
% axis equal; axis off;
% title('contours of estimated neurons');
% 
% plot contours with IDs
% [Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, T, min(T, 1000)))));
% figure;
% Cn = imresize(Cn, [d1, d2]); 
% plot_contours(neuron.A, Cn, 0.8, 0, [], neuron.Coor, 2);
% colormap winter;
% title('contours of estimated neurons');

%% check spatial and temporal components by playing movies
% save_avi = false;
% avi_name = 'play_movie.avi';
% neuron.Cn = Cn;
% neuron.runMovie(Ysignal, [0, 50], save_avi, avi_name);
% 
% %% save video
% kt = 3;     % play one frame in every kt frames
% save_avi = true;
% center_ac = median(max(neuron.A,[],1)'.*max(neuron.C,[],2)); % the denoised video are mapped to [0, 2*center_ac] of the colormap 
% cnmfe_save_video;

%% save results
globalVars = who('global');
eval(sprintf('save teststring %s',strjoin(globalVars)));
fprintf('Whole thing saved and done.');

