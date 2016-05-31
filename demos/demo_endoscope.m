%% clear workspace
clear; clc; close all;

%% choose data
if ~exist('nam', 'var') || isempty(nam)
    try
        load .dir.mat;
    catch
        dir_nm = [cd(), filesep];
    end
    [file_nm, dir_nm] = uigetfile(sprintf('%s*.tif', dir_nm));
    if dir_nm~=0
        save .dir.mat dir_nm;
    end
    nam = [dir_nm, file_nm];
end

temp = imread(nam);
[d1,d2, ~] = size(temp);
neuron_raw = Sources2D('d1',d1,'d2',d2);                         % dimensions of datasets
neuron_raw.Fs = 10;
neuron_raw.updateParams('ssub', 1, 'tsub', 3, ...
    'gSig', 4, 'gSiz', 15, 'bSiz', 1, ...
    'search_method', 'dilate', ...
    'merge_thr', 0.85, 'bas_nonneg', 1);    %% it seems that the
% raw data has been downsampled by ICA. adjust the downsampling factor
% here

%% load data
tic;
sframe=1;						% user input: first frame to read (optional, default 1)
num2read= 9000;					% user input: how many frames to read   (optional, default until the end)
[Y, neuron] = neuron_raw.load_data(nam, sframe, num2read);
[d1,d2, T] = size(Y);
Y = neuron.reshape(Y, 1);
neuron_raw.P.p = 2;      %order of AR model
fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);

%% initialization of A, C
tic;
debug_on = true;
save_avi = false;
K = 500; % maximum number of neurons to search. you can use [] to search the number automatically
neuron.options.min_corr = 0.7;
[Ain, Cin, center, Cn] = greedyROI_inscopix(Y, K, neuron.options, debug_on, save_avi);
fprintf('Time cost in computing the correlation image:     %.2f seconds\n', toc);
neuron.Cn = Cn;

figure;
imagesc(Cn); colorbar; axis equal off tight;
hold on;
plot(center(:, 2), center(:, 1), '*r');

[~, srt] = sort(max(Cin, [], 2).*max(Ain', [], 2), 'descend');
Ain = Ain(:, srt);
Cin = Cin(srt, :);

% the time cost in this step is highly dependent on the computer. It could
% be fast if you have a multicore CPU, the more the faster.
neuron.A = Ain;
neuron.C = Cin;

% %deconvole all temporal components
% C_med = quantile(Cin, 0.1,  2);
% Cin = bsxfun(@minus, Cin, C_med);
% neuron.C = Cin;

% neuron.deconvTemporal();
% neuron.displayNeurons([], Cin); %, 'garret_day1/neurons');

%% udpate background (block 1, the following three blocks can be run iteratively)
tic;

Ybg = Y-neuron.A*neuron.C;
nb = 10;     % rank of the background
neuron.updateBG(Ybg, nb, 'svd');  % here we use SVD to model the background
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
%% update spatial components (blcok 2)
tic;
max_min_ratio = 15;     % for each neuron's spatial component, it threshold the nonzero pixels to be bigger than max / max_min_ratio.
neuron.trimSpatial(max_min_ratio);
ind_nonzero = (neuron.A>0);     % nonzero pixels
neuron.options.se = strel('disk', 5);
IND = determine_search_location(neuron.A, 'dilate', neuron.options);

% update spatial components with model DY = A*DC
DY = diff(Y-neuron.b*neuron.f, 1, 2);
DC = diff(neuron.C, 1, 2);
DY(bsxfun(@lt, abs(DY), 2*std(DY, 0, 2))) = 0;
DC(bsxfun(@lt, abs(DC), 2*std(DC, 0, 2))) = 0;
tic; A = HALS_spatial(DY, neuron.A, DC, IND, 50);
% update spatial components with model Y = A*C;
% tic; A = HALS_spatial(Ysignal, neuron.A, neuron.C, IND, 100);
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

%% pick more neurons from the residual 
Yres = Ysignal - neuron.A * neuron.C; 
% neuron.auto_add(Yres); 
neuron.manual_add(Yres); 
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
