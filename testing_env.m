%% Test Script for optimizing CNMF_E parameters

%% clear workspace
clear; clc; close all;  
global  data d1 d2 numFrame ssub tsub sframe num2read Fs Ysiz Ybg_weights; % global variables, don't change them manualyl


%% choose data and map to memory

dir_nm = [cd(), filesep];
[file_nm, dir_nm] = uigetfile(fullfile(dir_nm, '*.tif;*.mat'));
nam = [dir_nm, file_nm];  % full name of the data file
[dir_nm, file_nm, file_type] = fileparts(nam);
nam_mat = [dir_nm, filesep, file_nm, '.mat'];
if strcmpi(file_type, '.mat')
    fprintf('The selected file is *.mat file\n');
elseif  exist(nam_mat, 'file')
    % the selected file has been converted to *.mat file already
    fprintf('The selected file has been replaced with its *.mat version.\n');
elseif or(strcmpi(file_type, '.tif'), strcmpi(file_type, '.tiff'))
    % convert
    tic;
    fprintf('converting the selected file to *.mat version...\n');
    nam_mat = tif2mat(nam);
    fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);
end

data = matfile(nam_mat);
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames

fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));

%% load in data from .mat file

sframe = 1; %starting frame
num2read = numFrame; %total number of frames to use
ssub = 1;
tsub = 1;

options.ssub = ssub;
options.tsub = tsub;

if and(ssub==1, tsub==1)
    Y = double(data.Y(:,:,sframe+(1:num2read)-1));
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else
    Tbatch = round((2^28)/(d1*d2)/tsub)*tsub;
    [~, ~, file_type2] = fileparts(nam_mat);
    temp_img = data.Y(:, :, 1);
    num2read = min(num2read, numFrame-sframe+1);
    if Tbatch>=num2read
        % load all data because the file is too small
        Yraw = data.Y(:, :, (1:num2read)+sframe-1);
        Y = dsData(double(Yraw),options);
    else
        % load data in batch model
        [d1s, d2s] = size(imresize(temp_img, 1/ssub));  % size of spatial downsampled data
        Ts = floor(num2read/tsub);  % frame number after downsampling
        Y = zeros(d1s, d2s, Ts);    % downsampled data size
        lastframe = sframe + Ts*tsub -1;  % index of the last frame for loading
        frame0 = sframe;
        while sframe<= lastframe
            tmp_num2read =  min(lastframe-sframe+1, Tbatch);
            fprintf('load data from frame %d to frame %d of %d total frames\n', sframe, sframe+tmp_num2read-1, lastframe);
            Yraw = data.Y(:, :, sframe:(sframe+tmp_num2read-1));
            temp = dsData(double(Yraw),options);
            Y(:, :, (sframe-frame0)/tsub + (1:size(temp, 3))) = temp;
            sframe = sframe + size(Yraw,3);
        end
    end
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n',...
        d1s, d2s, T, d1s*d2s*T*8/(2^30));
end

Y = reshape(Y,prod([d1 d2]),[]); % reshapes 3D video (tensor) into 2D matrix (pixels x time)

%% Create correlation image for initialization

options.d1 = d1;
options.d2 = d2;
options.gSig = 6;
options.gSiz = 15;
options.noise_range = [0.45 0.5];
options.noise_method = 'logmexp';
options.block_size = [64 64];

[Cn, PNR] = correlation_image_endoscope(Y, options);

%% the following parameters are required by the 'greedyROI_endoscope.m' initialization function
options.seed_method = 'auto'; %or manual
options.min_corr = 0.65;
options.min_pnr = 9;
options.min_pixel = 5;
options.deconv_flag = 1;
options.merge_thr = 0.5;
options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', false, ... % optimize the baseline
    'optimize_smin', true);  % optimize the threshold

tau_decay = 2;  %
tau_rise = 0.1;
nframe_decay = ceil(10*tau_decay*options.Fs);  % number of frames in decaying period
bound_pars = false;     % bound tau_decay/tau_rise or not
options.kernel = create_kernel('exp2', [tau_decay, tau_rise]*options.Fs, nframe_decay, [], [], bound_pars);


%% detrend (optional) and use greedy initialization to find cell peaks

nk = 2; %number of knots in bspline basis

[Ydt, ~, ~] = detrend_data(Y,nk);

[results, center, Cn, PNR, save_avi] = greedyROI_endoscope(Y, 200, options,true, false);


