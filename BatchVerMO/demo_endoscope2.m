function [A0s,File]=demo_endoscope2(bg_neuron_ratio,merge_thr,with_dendrites,K,start_frame,num_2read,name,neuron_full_partial,Mode,picname,File,Afinal,thresh_detecting_frames)
%function [A0s,File]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,min_pixel,bd,FS,SSub,TSub,bg_neuron_ratio,name,Mode,picname,Afinal,File,convolveType,merge_thr)
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
mode=Mode;                                 % 'initiation' mode or 'massive' mode
if strcmp(mode,'initiation')
    nam=name;
    cnmfe_choose_data;
elseif strcmp(mode,'massive')
    nam=name{1};
    display(nam)
    cnmfe_choose_data;
    Ysignal=name{2};
end
neuron_full=neuron_full_partial;
clear neuron_full_partial
neuron_full.updateParams('d1',d1, 'd2',d2);
min_pixel=neuron_full.options.min_pixel;
min_pnr=neuron_full.options.min_pnr;
min_corr=neuron_full.options.min_corr;

%% create Source2D class object for storing results and parameters
Fs = neuron_full.Fs;                       % frame rate
ssub = neuron_full.options.ssub;           % spatial downsampling factor
tsub = neuron_full.options.tsub;           % temporal downsampling factor

% with dendrites or not 
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron_full.options.search_method = 'dilate'; 
    neuron_full.options.bSiz = 20;
else
    % determine the search locations by selecting a round area
    neuron_full.options.search_method = 'ellipse';
    neuron_full.options.dist = 5;
end


%% downsample data for fast and better initialization
sframe=start_frame;					% user input: first frame to read (optional, default:1)
if isempty(num_2read)
    num2read = numFrame;             % user input: how many frames to read   (optional, default: until the end)
else
    num2read = num_2read;
end

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

neuron.options.nk = 1;  % number of knots for detrending 

% greedy method for initialization
tic;
if strcmp(mode,'initiation')
    display(['working on ',nam])
    [center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
    fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);
    if isempty(neuron.A)
        A0s=neuron.A;
        File.options=[]; File.Ysignal=[]; File.neuron=[]; File.Ybg=[];
        clear global
        return
    end
    % sort neurons
    [~, srt] = sort(max(neuron.C, [], 2), 'descend');
    neuron.orderROIs(srt);
    neuron_init = neuron.copy();
        
elseif strcmp(mode,'massive')
    % parameters, estimate the background
    spatial_ds_factor = 1;              % spatial downsampling factor. it's for faster estimation
    thresh = 10;                        % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
    [C,~]=extract_c(Ysignal,[],Afinal);
    neuron.A=Afinal;
    neuron.C_raw=C;
    neuron.C=C;
    [~,~]=neuron.updateTemporal_endoscope(Ysignal,false);
    cnmfe_update_BG;
    [~,~]=neuron.updateTemporal_endoscope(Ysignal,false);
    A0s=[];
    File.neuron=neuron.copy();
    clear global
    return 
end

%% iteratively update A, C and B
% parameters, merge neurons
display_merge = false;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the background
spatial_ds_factor = 1;          % spatial downsampling factor. it's for faster estimation
thresh=thresh_detecting_frames; % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
%bg_neuron_ratio = bg_neuron_ratio;  % spatial range / diameter of neurons

% parameters, estimate the spatial components
update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
Nspatial = 5;       % this variable has different meanings: 
                    %1) udpate_spatial_method=='hals' or 'hals_thresh',
                    %then Nspatial is the maximum iteration 
                    %2) update_spatial_method== 'nnls', it is the maximum
                    %number of neurons overlapping at one pixel 
               
% parameters for running iterations 
nC = size(neuron.C, 1);    % number of neurons 

maxIter = 2;        % maximum number of iterations 
miter = 1; 
while miter <= maxIter
    if strcmp(mode,'initiation')
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
        if isempty(neuron.A); A0s=neuron.A; File.options=[]; File.Ysignal=[]; File.neuron=[]; File.Ybg=[];
            clear global; return; end
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
        if strcmp(mode,'initiation')
            neuron.updateTemporal_endoscope(Ysignal,true);
            if isempty(neuron.A); break; end
            cnmfe_quick_merge;              % run neuron merges
            if isempty(neuron.A); break; end
        elseif strcmp(mode,'massive')
            neuron.updateTemporal_endoscope(Ysignal,false);
        end
    %spatial
        if strcmp(mode,'initiation')
            neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method,[],true);
            if isempty(neuron.A); break; end
            neuron.trimSpatial(0.01, 3, min_pixel); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
            if isempty(neuron.A); break; end
            if isempty(merged_ROI)
                break;
            end
        elseif strcmp(mode,'massive')
            neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method,[],false);
            neuron.trimSpatial(0.01, 3, min_pixel, false);
        end
    end
    fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);    
    
    %% pick neurons from the residual (cell 4).
    if strcmp(mode,'initiation')
        if miter==1
            neuron.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
            [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par); % method can be either 'auto' or 'manual'
        end
        if isempty(neuron.A); A0s=neuron.A; File.options=[]; File.Ysignal=[]; File.neuron=[]; File.Ybg=[];
            clear global; return; end
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

%if isempty(neuron.A); A0s=neuron.A; File.options=[]; File.Ysignal=[]; clear global; return; end
%Ybg=Ybg+b0;
%Ysignal_sn=Ysignal;
%noise=neuron.P.sn_neuron;
% if strcmp(mode,'initiation'); Ysignal=neuron.A*neuron.C_raw; end

%% apply results to the full resolution
if or(ssub>1, tsub>1)
    neuron_ds = neuron.copy();  % save the result
    neuron = neuron_full.copy();
    cnmfe_full;
    neuron_full = neuron.copy();
end

%% for 'massive' mode,  save results
if strcmp(mode,'massive')
    A0s=[];
    File.neuron=neuron;
end
%% for 'initiation' mode see and save results
resultstring=sprintf('%s_results', Picname);
neuron.viewNeurons([], neuron.C_raw, resultstring);
close(gcf);

neuron.drawPNRCn(min_pnr,min_corr)
close(gcf);

ColorAllNeurons(neuron.A,d1,d2,[Picname,strcat('PNR=',num2str(min_pnr)),'.png'],outputdir);
if strcmp(mode,'initiation')
    A0s=neuron.A;
    File.options=neuron.options;    
    File.neuron=neuron.copy();
    File.Ysignal=Ysignal;
    File.Ybg=Ybg;
end
clear global
%globalVars = who('global');
%eval(sprintf('save %s%s%s_results.mat %s', dir_nm, filesep, file_nm, strjoin(globalVars)));