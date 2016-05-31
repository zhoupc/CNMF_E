%% clear workspace
clear; clc; close all;

%% select data and map it with the memory
if ~exist('nam', 'var') || isempty(nam)
    try
        load .dir.mat; %load previous path
    catch
        dir_nm = [cd(), filesep]; %use the current path
    end
    [file_nm, dir_nm] = uigetfile(fullfile(dir_nm, '*.tif;*.mat'));
    if dir_nm~=0
        save .dir.mat dir_nm;
    else
        fprintf('no file was selecteD. STOP!\N');
        return;
    end
    nam = [dir_nm, file_nm];  % full name of the data file
    [~, file_nm, file_type] = fileparts(nam);
end

% convert the data to mat file
nam_mat = [dir_nm, file_nm, '.mat']; 
if strcmpi(file_type, '.mat')
    fprintf('The selected file is *.mat file\n'); 
elseif  exist(nam_mat', 'file')
    % the selected file has been converted to *.mat file already 
    fprintf('The selected file has been replaced with its *.mat version\n'); 
elseif or(strcmpi(file_type, '.tif'), strcmpi(file_type, '.tiff'))
    % convert 
    tic;
    fprintf('converting the selected file to *.mat version...\n'); 
    nam_mat = tif2mat(nam);
    fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);
else
    fprintf('The selected file type was not supported yet! email me to get support (zhoupc1988@gmail.com)\n');
    return; 
end

data = matfile(nam_mat); 
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames

fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));

%% create Source2D class object for storing results and parameters
neuron_raw = Sources2D('d1',d1,'d2',d2);   % dimensions of datasets
neuron_raw.Fs = 10;         % frame rate
neuron_raw.options.Fs = 10; 
neuron_raw.options.trend_itval = 100;   % interval between knots of the spline basis 
ssub = 1;           % spatial downsampling factor 
tsub = 1;           % temporal downsampling factor 
neuron_raw.updateParams('ssub', ssub,...  % spatial downsampling factor
    'tsub', tsub, ...  %temporal downsampling factor
    'gSig', 4,... %width of the gaussian kernel, which can approximates the average neuron shape
    'gSiz', 15, ...% average size of a neuron
    'bSiz', 2, ... % half width of the kernel to dilate nonzero pixels of each neuron
    'search_method', 'dilate', ... % searching method
    'merge_thr', 0.7, ... % threshold for merging neurons
    'bas_nonneg', 1);   % 1: positive baseline of each calcium traces; 0: any baseline

%% downsample data for fast and better initialization 
sframe=1;						% user input: first frame to read (optional, default:1)
num2read= numFrame;             % user input: how many frames to read   (optional, default: until the end)

tic; 
if and(ssub==1, tsub==1)
    neuron = neuron_raw; 
    Y = data.Y;  
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else 
    [Y, neuron] = neuron_raw.load_data(nam_mat, sframe, num2read);
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
end 
Y = neuron.reshape(Y, 1);
neuron_raw.P.p = 2;      %order of AR model

fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);

%% initialization of A, C
tic;
debug_on = true;
save_avi = true;
neuron.options.min_corr = 0.7;
neuron.options.nk = 0; %round(T/(180*neuron.Fs)); % number of knots for spline basis, the interval between knots is 100 seconds
patch_par = 1; %[2,2];  % divide the optical field into 4 X 4 patches and do initialization patch by patch
K = 500; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
neuron.options.bd = 20; 
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); 
neuron_init = neuron.copy(); 
neuron.options.merge_thr = 0.65;     % threshold for merging neurons
merged_ROI = neuron.quickMerge('C');  %merge neurons based on their temporal correlation 

% display merged neurons 
figure; 
for m=1:length(merged_ROI)
    subplot(221); 
    neuron.image(sum(neuron_init.A(:, merged_ROI{m}), 2)); 
    axis equal off tight; 
    subplot(2,2,3:4); 
    plot(neuron_init.C(merged_ROI{m}, :)'); 
    axis tight; 
    pause; 
end

%% order neurons and delete some low quality neurons 
A = neuron.A; 
C = neuron.C; 
% [~, srt] = sort(max(A, [], 1)'.*max(C, [], 2), 'descend'); 
[Cpnr, srt] = sort(max(C, [], 2)./get_noise_fft(C), 'descend'); 
neuron.orderROIs(srt); 
A = neuron.A; 
C = neuron.C; 
neuron.viewNeurons(); 

% display contours of the neurons 
figure; 
neuron.viewContours(Cn, 0.9, 0); 

% deconvolve all traces 
% Cin = neuron.C; 
% neuron.deconvTemporal();
% neuron.displayNeurons([], Cin);

%% udpate background (block 1, the following three blocks can be run iteratively)
tic;

Ybg = double(Y)-neuron.A*neuron.C;
nb = 10;     % rank of the background
neuron.updateBG(Ybg, nb, 'svd');  % here we use SVD to model the backgroun)
clear Ybg;
fprintf('Time cost in inferring the background:     %.2f seconds\n', toc);

Ysignal = Y-neuron.b*neuron.f; % data after removing the background
% neuron.playMovie(Ysignal); % play the movie after subtracting the background.
figure('position', [100, 100, 1000, 350]);
for m=1:nb
    subplot(131);
    neuron.image(neuron.b(:, m));
    axis equal off tight;
    subplot(1,3,2:3);
    plot(neuron.f(m, :));
    pause;
end

%% pick neurons from the residual 
neuron.manual_add(Ysignal-neuron.A*neuron.C); 
%% update spatial components (blcok 2)
tic;
max_min_ratio = 15;     % for each neuron's spatial component, it threshold the nonzero pixels to be bigger than max / max_min_ratio.
neuron.trimSpatial(max_min_ratio);
ind_nonzero = (neuron.A>0);     % nonzero pixels
neuron.options.se = strel('disk', 5);
IND = determine_search_location(neuron.A, 'dilate', neuron.options);

% update spatial components with model DY = A*DC
% DY = diff(Y-neuron.b*neuron.f, 1, 2);
% DC = diff(neuron.C, 1, 2);
% DY(bsxfun(@lt, abs(DY), 2*std(DY, 0, 2))) = 0;
% DC(bsxfun(@lt, abs(DC), 2*std(DC, 0, 2))) = 0;
% tic; A = HALS_spatial(DY, neuron.A, DC, IND, 50);
% update spatial components with model Y = A*C;
tic; A = HALS_spatial(Ysignal, neuron.A, neuron.C, IND, 100);
A = full(A);
ind_del = false(1, size(A, 2));
for m=1:size(A,2)
    tmp = neuron.reshape(A(:,m), 2);
    l = bwlabel(tmp>(max(A(:, m))/max_min_ratio), 4);
    label_seed = mode(l(ind_nonzero(:, m)));
    if label_seed==0
        ind_del(m) = true;
    else
        tmp(l~=label_seed) = 0;
    end
    A(:, m) = tmp(:);
end
neuron.A = A;
neuron.delete(ind_del);
fprintf('Time cost in updating neuronal spatial components:     %.2f seconds\n', toc);

%% update C  (block 3)
tic;
% A = neuron.A;
% C = (A'*A)\(A'*Ysignal);

C = HALS_temporal(Ysignal, neuron.A, neuron.C, 100);
neuron.C = C;
% neuron.deconvTemporal();    %deconvolve all temporal component if you want

fprintf('Time cost in updating neuronal temporal components:     %.2f seconds\n', toc);

%% display neurons
figure;
neuron.viewNeurons();

%% view contours
figure;
neuron.viewContours(Cn, 0.8, 0);
axis equal; axis off;
title('contours of estimated neurons');

%% check results by visualizing all spatial and temporal components one by one
folder_nm = [];%'neurons';
neuron.viewNeurons([], C, folder_nm);

%% check spatial and temporal components by playing movies
save_avi = true;
avi_name = 'play_movie.avi';
neuron.runMovie(Ysignal, [0, 100], save_avi, avi_name);
%%
