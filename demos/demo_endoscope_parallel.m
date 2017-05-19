% Demo Script for running CNMF_E in parallel
% Based on Pengcheng Zhou's CNMF_E demo script, 'demo_endoscope,' with
% simple parallelization implemented.
% Conor Heins, conor.heins@gmail.com
% May 18 2017

%% clear workspace
clear; clc; close all;  
global  d1 d2 numFrame ssub tsub sframe num2read Fs neuron neuron_ds ...
    neuron_full; %#ok<NUSED> % global variables, don't change them manually

%% select data and map it to RAM
% nam = 'E:\CNMF_E-master\demos\data_endoscope.tif';
cnmfe_choose_data;

%% create Source2D class object for storing results and parameters
Fs = 6;             % frame rate
ssub = 2;           % spatial downsampling factor
tsub = 2;           % temporal downsampling factor
gSig = 3;           % width of the gaussian kernel, which can approximates the average neuron shape
gSiz = 10;          % maximum diameter of neurons in the image plane. larger values are preferred.
neuron_full = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'ssub', ssub, 'tsub', tsub, ...  % downsampleing
    'gSig', gSig,...    % sigma of the 2D gaussian that approximates cell bodies
    'gSiz', gSiz,...    % average neuron size (diameter)
    'use_parallel',false,... %disable parallellization within CNMF_E to avoid transparency violations
    'temporal_parallel',false); %disable parallellization within CNMF_E to avoid transparency violations     
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

%% load small portion of data for displaying correlation image
% balance number of frames to read against the image dimensions (d1 x d2) to
% limit RAM footprint
sframe = 1;
num2read = 1000; 

if and(ssub==1, tsub==1)
    neuron_small = neuron_full;
    Y = double(data.Y(:, :, sframe+(1:num2read)-1));
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else
    [Y, neuron_ds] = neuron_full.load_data(nam_mat, sframe, num2read);
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
    neuron_small = neuron_ds.copy(); 
end

Y = neuron_small.reshape(Y,1);

%% compute correlation image and peak-to-noise ratio for selected portion of data
% this step is not necessary, but it can give you some hints on parameter selection, e.g., min_corr & min_pnr

[Cn, pnr] = neuron_small.correlation_pnr(Y(:, round(linspace(1, T, min(T, 1000)))));

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
%% create indices for splitting field-of-view into spatially-overlapping patches (for parallel processing)

patch_size = [64, 64]; %patch size
overlap = [12 12]; %patch overlap
min_patch_sz = [16 16]; %minimum patch size in either direction
patches = construct_patches(Ysiz(1:end-1),patch_size,overlap,min_patch_sz);

%% Load and run CNMF_E on full dataset in patches

sframe=1;						% user input: first frame to read (optional, default:1)
num2read= numFrame;             % user input: how many frames to read (optional, default: until the end)

RESULTS(length(patches)) = struct();

%%  PARALLEL CNMF_E
% If you change main loop to a for loop (sequential processing), you can save space by condensing many of the 
% steps below into scripts, as is done in the original demo_endoscope.m

parfor i = 1:length(patches) 
    
    % Load data from individual (i-th) patch and store in temporary Sources2D() object ('neuron_patch')                      
                             
    neuron_patch = neuron_full.copy();
    sframe_patch = sframe; num2read_patch = num2read;
    
    if and(ssub==1, tsub==1)
        Y = double(data.Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),sframe_patch+(1:num2read_patch)-1));
        [d1p, d2p, T] = size(Y);
    else
        Yraw = data.Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),sframe_patch+(1:num2read_patch)-1);
        [neuron_patch_ds, Y] = neuron_patch.downSample(double(Yraw));
        [d1p, d2p, T] = size(Y);
        neuron_patch = neuron_patch_ds.copy(); 
    end

    neuron_patch.options.d1 = d1p; neuron_patch.options.d2 = d2p;
    Y = neuron_patch.reshape(Y, 1);       % convert a 3D video into a 2D matrix
    
    %% Initialization of A,C -- parameters
    % must initialize individual update-parameters within parfor loop to avoid
    % transparency violations
    
    debug_on = false;   % visualize the initialization procedue.
    save_avi = false;   %save the initialization procedure as an avi movie. 
    patch_par = [1,1]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large 
    K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
    
    min_corr = 0.85;     % minimum local correlation for a seeding pixel
    min_pnr = 20;       % minimum peak-to-noise ratio for a seeding pixel
    min_pixel = 5;      % minimum number of nonzero pixels for each neuron
    bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
    neuron_patch.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
        'min_pixel', min_pixel, 'bd', bd);
    neuron_patch.options.nk = 1;  % number of knots for detrending 
    
    % greedy method for initialization
    neuron_patch.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
    
    %% iteratively update A, C and B
    % must initialize individual update-parameters within parfor loop to avoid
    % transparency violations
    
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
    
    % parameters for running iterations
    nC = size(neuron_patch.C, 1);    % number of neurons
    
    
    %% Continue with full CNMF_E if given patch has any neurons
    
    if nC>0 % check for presence of neurons in patch
        % sort neurons
        [~, srt] = sort(max(neuron_patch.C, [], 2), 'descend');
        neuron_patch.orderROIs(srt);
        
        maxIter = 2;        % maximum number of iterations 
        miter = 1;
        while miter <= maxIter
           %% merge neurons, order neurons and delete some low quality neurons
           
           if nC <= 1 % check if there's at most 1 neuron
               break
           end
           
           if miter ==1
               merge_thr = [1e-1, 0.8, .1];     % thresholds for merging neurons
               % corresponding to {sptial overlaps, temporal correlation of C,
               %temporal correlation of S}
           else
               merge_thr = [0.6, 0.5, 0.1];
           end

           % merge neurons
           neuron_patch.quickMerge(merge_thr); % run neuron merges
           %sort neurons
           [~,srt] = sort(max(neuron_patch.C,[],2).*max(neuron_patch.A,[],1)','descend');
           neuron_patch.orderROIs(srt);
            
           %% udpate background
           tic;
           
           Ybg = Y-neuron_patch.A*neuron_patch.C;
           rr = ceil(neuron_patch.options.gSiz * bg_neuron_ratio);
           active_px = []; %(sum(IND, 2)>0); %If some missing neurons are not covered by active_px, use [] to replace IND
           [Ybg, ~] = neuron_patch.localBG(Ybg,spatial_ds_factor,rr,active_px,neuron_patch.P.sn,thresh); %estimate local background
           
           %subtract background from the raw data to obtain signal for
           %subsequent CNMF
           Ysignal = Y - Ybg;
           
           % estimate noise
           if ~isfield(neuron_patch.P,'sn') || isempty(neuron_patch.P.sn)
               % estimate the noise for all pixels
               b0 = zeros(size(Ysignal,1), 1);
               sn = b0;
               for m = 1:size(neuron_patch.A,1)
                   [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m,:));
               end
               Ysigma = bsxfun(@minus,Ysignal,b0);
           end
           
           fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);           
            %% update spatial & temporal components
            tic;
            for m=1:2
                %temporal
                neuron_patch.updateTemporal_endoscope(Ysignal);
                
                % merge neurons
                [merged_ROI, ~] = neuron_patch.quickMerge(merge_thr); % run neuron merges
                %sort neurons
                [~,srt] = sort(max(neuron_patch.C,[],2).*max(neuron_patch.A,[],1)','descend');
                neuron_patch.orderROIs(srt);
                
                %spatial
                neuron_patch.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
                neuron_patch.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
                if isempty(merged_ROI)
                    break;
                end
            end
            fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
            
            %% pick neurons from the residual (cell 4).
            if miter==1
                neuron_patch.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
                [center_new, Cn_res, pnr_res] = neuron_patch.pickNeurons(Ysignal - neuron_patch.A*neuron_patch.C, patch_par); % method can be either 'auto' or 'manual'
            end
            
            %% stop the iteration
            temp = size(neuron_patch.C, 1);
            if or(nC==temp, miter==maxIter)
                break;
            else
                miter = miter+1;
                nC = temp;
            end
        end
    end
    
    %% Store results from individual patch in master structure with every patch's output
    
    RESULTS(i).A = neuron_patch.A;
    RESULTS(i).C = neuron_patch.C;
    RESULTS(i).C_raw = neuron_patch.C_raw;
    RESULTS(i).S = neuron_patch.S;
    RESULTS(i).P = neuron_patch.P;
    fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(length(patches)), '.\n']);
    
end

%% Store all patches in full Sources2D object for full dataset 

neuron_full.P.kernel_pars = [];
for i = 1:length(patches)
    if size(RESULTS(i).A,2)>0
        Atemp = zeros(neuron_full.options.d1,neuron_full.options.d2,size(RESULTS(i).A,2));
        for k = 1:size(RESULTS(i).A,2)
            Atemp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),k)=reshape(RESULTS(i).A(:,k),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
        end
        neuron_full.A = [neuron_full.A,reshape(Atemp,d1*d2,k)];
        neuron_full.C = [neuron_full.C;RESULTS(i).C];
        neuron_full.C_raw = [neuron_full.C_raw;RESULTS(i).C_raw];
        neuron_full.S = [neuron_full.S;RESULTS(i).S];
        neuron_full.P.sn = [neuron_full.P.sn,RESULTS(i).P.sn];
        neuron_full.P.kernel_pars = [neuron_full.P.kernel_pars;RESULTS(i).P.kernel_pars];
    end
    clear Atemp;
end

neuron = neuron_full.copy();

%% delete, trim, split neurons
neuron.viewNeurons([], neuron.C_raw); 

%% merge neurons
display_merge = true; % set to true if you want to inspect every candidate merge
view_neurons = false; % set to true if you want to inspect all neurons after quick merge routine

merge_thr = [0.6, 0.5, 0.1];  % choose thresholds for merging neurons (this will primarily merge neurons redundantly found by multiple patch processes, likely in the patch-overlaps)
cnmfe_quick_merge;            % run neuron merges

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
neuron.Coor = neuron.get_contours(0.8); % energy within the contour is 80% of the total 

% plot contours with IDs
figure;
Cn = imresize(Cn, [d1, d2]); 
plot_contours(neuron.A, Cn, 0.8, 0, [], neuron.Coor, 2);
colormap winter;
title('contours of estimated neurons');
neuron.Cn = Cn;

%% save results
globalVars = who('global');
eval(sprintf('save %s%s%s_results.mat %s', dir_nm, filesep, file_nm, strjoin(globalVars)));