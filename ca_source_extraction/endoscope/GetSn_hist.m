function sn = GetSn_hist(Y, use_parallel)
%% Estimate noise level for each row. It fits the histogram of row values using a gaussian function

%% inputs:
%   Y: N X T matrix, fluorescence trace

%% outputs:
%   sn: scalar, std of the noise

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016

%% input arguments
if any(size(Y)==1)
    Y = reshape(Y, 1, []);
end
[nrows, ncols] = size(Y);
if ~exist('use_parallel', 'var') || isempty(use_parallel)
    use_parallel = true;
end
%% estimate the noise levels for each pixel
sn = cell(nrows, 1);
Y = mat2cell(Y, ones(1, nrows), ncols);
if use_parallel
    parfor m=1:nrows
        [~, sn{m}] = estimate_baseline_noise(Y{m});
    end
else
    for m=1:nrows
        [~, sn{m}] = estimate_baseline_noise(Y{m});
    end
end
sn = cell2mat(sn);
