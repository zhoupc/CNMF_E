function [Cn, PNR] = correlation_pnr_parallel(obj, frame_range)
%% compute correlation image and pnr imaging in parallel
%% input:
%   frame_range: 1 X 2 vector indicating the starting and ending frames

%% Output:
%   Cn:     correlation image
%   PNR:    peak to noise ratio
%% Author: Pengcheng Zhou, Columbia University, 2017
%% email: zhoupc1988@gmail.com

%% process parameters

% map data
mat_data = obj.P.mat_data;

% dimension of data
dims = mat_data.dims;
d1 = dims(1);
d2 = dims(2);
T = dims(3);
obj.options.d1 = d1;
obj.options.d2 = d2;

% parameters for patching information
patch_pos = mat_data.patch_pos;
block_pos = mat_data.block_pos;

% number of patches
[nr_patch, nc_patch] = size(patch_pos);

% downsampling
ssub = obj.options.ssub;
tsub = obj.options.tsub;

% frames to be loaded for initialization
if ~exist('frame_range', 'var')
    frame_range = obj.frame_range;
end
if isempty(frame_range)
    frame_range = [1, T];
else
    frame_range(frame_range<1) = 1;
    frame_range(frame_range>T) = T;
end
T = diff(frame_range) + 1;
obj.frame_range = frame_range;

% parameters for detrending
if isfield(obj.options, 'nk') % number of knots for creating spline basis
    nk = obj.options.nk;
else
    nk = 1;
end
detrend_method = obj.options.detrend_method;

% parameter for avoiding using boundaries pixels as seed pixels
options = obj.options;

%% compute correlation image and pnr image in parallel
Cn = zeros(d1, d2);
PNR = zeros(d1, d2);

results = cell(nr_patch, nc_patch);

for m=1:(nr_patch*nr_patch)
    fprintf('|');
end
fprintf('\n');
parfor mpatch=1:(nr_patch*nc_patch)
    % get the indices corresponding to the selected patch
    tmp_patch = patch_pos{mpatch};
    tmp_block = block_pos{mpatch};
    
    % boundaries pixels to be avoided for detecting seed pixels
    tmp_options = options;
    
    % patch dimension
    tmp_d1 = diff(tmp_block(1:2))+1;
    tmp_d2 = diff(tmp_block(3:4))+1;
    tmp_options.gSig = tmp_options.gSig/ssub; 
    tmp_options.gSiz = ceil(tmp_options.gSiz/ssub); 
    
    % load the patch data
    Ypatch = get_patch_data(mat_data, tmp_patch, frame_range, true);
    if (ssub~=1)||(tsub~=1)
        Ypatch = dsData(Ypatch, tmp_options);
    end
    [tmp_options.d1, tmp_options.d2, T] = size(Ypatch);
    Ypatch = double(reshape(Ypatch, [], T));
    if nk>1
        Ypatch_dt = detrend_data(Ypatch, nk, detrend_method); % detrend data
        [tmp_Cn, tmp_PNR] = correlation_image_endoscope(Ypatch_dt, tmp_options);
    else
        [tmp_Cn, tmp_PNR] = correlation_image_endoscope(Ypatch, tmp_options);
    end
    if (ssub~=1)
        tmp_Cn = imresize(tmp_Cn, [tmp_d1, tmp_d2]);
        tmp_PNR = imresize(tmp_PNR, [tmp_d1, tmp_d2]);
    end
    % put everthing into one struct variable
    results{mpatch} = struct('Cn', tmp_Cn, 'PNR', tmp_PNR);
    fprintf('.');
end
fprintf('\n');

%% collect results
for mpatch=1:(nr_patch*nc_patch)
    % get the indices corresponding to the selected patch
    [mr, mc] = ind2sub([nr_patch, nc_patch], mpatch);
    tmp_patch = patch_pos{mr, mc};
    tmp_block = block_pos{mr, mc};
    r0 = tmp_block(1);
    r1 = tmp_block(2);
    c0 = tmp_block(3);
    c1 = tmp_block(4);
    ind_patch = true(r1-r0+1, c1-c0+1);
    ind_patch((tmp_patch(1):tmp_patch(2))-r0+1, (tmp_patch(3):tmp_patch(4))-c0+1) = false;
    
    % unpack results
    tmp_results = results{mpatch};
    %     eval(sprintf('tmp_results=results_patch_%d;', mpatch));
    tmp_Cn = tmp_results.Cn;
    tmp_PNR = tmp_results.PNR;
    
    Cn(r0:r1, c0:c1) = max(Cn(r0:r1, c0:c1), tmp_Cn.*(1-ind_patch));
    PNR(r0:r1, c0:c1) = max(PNR(r0:r1, c0:c1), tmp_PNR.*(1-ind_patch));
end