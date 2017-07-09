function [A0s,File]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,min_pixel,bd,FS,SSub,TSub,bg_neuron_ratio,name,Mode,picname,Afinal,File,convolveType,merge_thr)
% [A0s,File]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,FS,SSub,TSub,bg_neuron_ratio,name,Mode,picname,Afinal,File)
% It is converted to function to cater to parfor loop.
%   Meanwhile, some part of cnmf-e demo is applied to this sampling process while some simplied methods are
%   applied for later for each file. Use 'mode' variable to specify this.
% In "initiation" mode, used in sampling data, most processes are the same as the original script,
% except that (1) initComponents_endoscope calls modified
%               greedy_ROI_endoscope() which plots PNR and corr for seeds
%               and non-neuron pixel. This will help you decide which PNR
%               to use.
%             (2) output is a File structure. File.options=neuron.options; File.Y=Y(raw signal); File.Ysignal=Ysignal; (Bg-subtracted and denoised Ysignal)
%                           and A0s, A0s is the output structure containing all A's obtained from each of the sample files.
% In "massive" mode, which is the process of applying the same A to all
% data in the folder, many steps are simplified.
%             (1) No initiation. Just get C from A using background
%             subtracted data.
%             (2) Then Deconvolve C.
%             In this mode, A0s is [], the File output only contains File.A=neuron.A; File.C=neuron.C; File.ind_del=ind_del;

% modified by Shijie Gu from demo_endoscope script by Pengcheng Zhou.
%% set global variables
global  d1 d2 numFrame ssub tsub sframe num2read Fs neuron neuron_ds ...
    neuron_full Ybg_weights mode Picname nam outputdir; %#ok<NUSED> % global variables, don't change them manually
Picname=picname;
%% select data and map it to the RAM
nam=name;
cnmfe_choose_data;

%% create Source2D class object for storing results and parameters

mode=Mode;          % 'initiation' mode or 'massive' mode
Fs = FS;             % frame rate
ssub = SSub;           % spatial downsampling factor
tsub = TSub;           % temporal downsampling factor
gSig = gSig;           % width of the gaussian kernel, which can approximates the average neuron shape
gSiz = gSiz;          % maximum diameter of neurons in the image plane. larger values are preferred.
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
neuron_full.options.deconv_options = struct('type', convolveType, ... % model of the calcium traces. {'ar1', 'ar2'}
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
% cnmfe_show_corr_pnr;    % this step is not necessary, but it can give you some...
                        % hints on parameter selection, e.g., min_corr & min_pnr

%% initialization of A, C
% parameters
debug_on = false;   % visualize the initialization procedue. 
save_avi = false;   %save the initialization procedure as an avi movie. 
patch_par = [1,1]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large 
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically

min_corr = min_corr;     % minimum local correlation for a seeding pixel
min_pnr = min_pnr;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = min_pixel;      % minimum number of nonzero pixels for each neuron
bd = bd;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);
neuron.options.nk = 1;  % number of knots for detrending 

% greedy method for initialization
tic;
if strcmp(mode,'initiation')
    display(['working on ',nam])
    [center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
    fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);
    if isempty(neuron.A)
        A0s=neuron.A;
        File.options=neuron.options;
        File.Ysignal=[];
        clear global
        return
    end
        
elseif strcmp(mode,'massive')
    % parameters, estimate the background
    spatial_ds_factor = 1;              % spatial downsampling factor. it's for faster estimation
    thresh = 10;                        % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
    bg_neuron_ratio = bg_neuron_ratio;  % spatial range / diameter of neurons
    BackgroundSub
    [C,~]=extract_c(Ysignal,[],Afinal);
    neuron.A=Afinal;
    neuron.C=C;
    [~,ind_del]=neuron.updateTemporal_endoscope(Ysignal,false);
    A0s=[];
%     File.A=neuron.A;
%     File.C=neuron.C;
    File.ind_del=ind_del;
    File.neuron=neuron;    
    return
elseif strcmp(mode,'BackgroundSubOnly')
    % parameters, estimate the background
    spatial_ds_factor = 1;              % spatial downsampling factor. it's for faster estimation
    thresh = 10;                        % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
    bg_neuron_ratio = bg_neuron_ratio;  % spatial range / diameter of neurons
    BackgroundSub
    return
end

% sort neurons
[~, srt] = sort(max(neuron.C, [], 2), 'descend');
neuron.orderROIs(srt);
neuron_init = neuron.copy();

%% iteratively update A, C and B
% parameters, merge neurons
display_merge = false;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the background
spatial_ds_factor = 1;      % spatial downsampling factor. it's for faster estimation
thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = bg_neuron_ratio;  % spatial range / diameter of neurons

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
while miter <= maxIter
    %% merge neurons, order neurons and delete some low quality neurons
     if miter ==1
        merge_thr = [1e-1, 0.8, .1];     % thresholds for merging neurons
        % corresponding to {sptial overlaps, temporal correlation of C,
        %temporal correlation of S}
    else
        merge_thr = merge_thr; 
    end
    % merge neurons
    cnmfe_quick_merge;              % run neuron merges
    if isempty(neuron.A)
        A0s=neuron.A;
        File.options=neuron.options;
        File.Ysignal=[];
        clear global
        return
    end
    
    %% udpate background (cell 1, the following three blocks can be run iteratively)
    % estimate the background
    tic;
    cnmfe_update_BG;
    fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
 
    %% update spatial & temporal components
    tic;
    for m=1:2  
        %temporal
        neuron.updateTemporal_endoscope(Ysignal,true);
        cnmfe_quick_merge;              % run neuron merges
        if isempty(neuron.A)
            A0s=neuron.A;
            File.options=neuron.options;
            File.Ysignal=[];
            clear global
            return
        end
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
%Ybg=Ybg+b0;
%Ysignal_sn=Ysignal;
%noise=neuron.P.sn_neuron;
Ysignal=neuron.A*neuron.C;

%% apply results to the full resolution
if or(ssub>1, tsub>1)
    neuron_ds = neuron.copy();  % save the result
    neuron = neuron_full.copy();
    cnmfe_full;
    neuron_full = neuron.copy();
end


%% delete some neurons and run CNMF-E iteration 
% neuron.viewNeurons([], neuron.C_raw); 
% tic;
% cnmfe_update_BG;
% fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
% %update spatial & temporal components
% tic;
% for m=1:2
%     %temporal
%     neuron.updateTemporal_endoscope(Ysignal);
%     cnmfe_quick_merge;              % run neuron merges
%     %spatial
%     neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
%     neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
% end
% fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

%% for 'initiation' mode see and save results
resultstring=sprintf('%s_results', Picname);
neuron.viewNeurons([], neuron.C_raw, resultstring);
close(gcf);

neuron.drawPNRCn(min_pnr,min_corr)
close(gcf);

ColorAllNeurons(neuron.A,d1,d2,Picname,outputdir);
if strcmp(mode,'initiation')
    A0s=neuron.A;
    File.options=neuron.options;
    File.Ysignal=Ysignal;
end
clear global
%globalVars = who('global');
%eval(sprintf('save %s%s%s_results.mat %s', dir_nm, filesep, file_nm, strjoin(globalVars)));