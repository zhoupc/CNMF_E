function [ai, ci, ind_success, sn] = extract_ac(HY, Y, ind_ctr, sz)
%% given a patch of raw & high-pass filtered calcium imaging data, extract
% spatial and temporal component of one neuron (ai, ci). if succeed, then
% return an indicator ind_succes with value 1; otherwise, 0.
%% inputs:
%       HY:     d X T matrix, filtered patch data
%       Y:      d X T matrix, raw data
%       ind_ctr:        scalar, location of the center
%       sz:         2 X 1 vector, size of the patch

%% Author: Pengcheng Zhou, Carnegie Mellon University.

%% parameters 
nr = sz(1);
nc = sz(2);
min_corr = 0.7;
min_pixels = 5;

%% find pixels highly correlated with the center
% HY(HY<0) = 0;       % remove some negative signals from nearby neurons
y0 = HY(ind_ctr, :);
tmp_corr = reshape(corr(y0', HY'), nr, nc);
data = HY(tmp_corr>min_corr, :);


%% estimate ci with the mean or rank-1 NMF
ci = mean(data, 1);
[~, sn] = estimate_baseline_sn(ci);  % estimate the noise level 
% ci = ci - min(ci); % avoid nonnegative baseline
% [~, ci] = nnmf(ci, 1);
if norm(ci)==0
    ai=[];
    ind_success=false;
    return;
end

%% extract spatial component
% estiamte the background level using the boundary
y_bg = median(Y(tmp_corr(:)<min_corr, :), 1); % using the median of the whole field (except the center area) as background estimation

% sort the data, take the differencing and estiamte ai
thr_noise = 5;      % threshold the nonzero pixels to remove noise
[~, ind_sort] = sort(y_bg, 'ascend');   % order frames to make sure the background levels are close within nearby frames
dY = diff(Y(:, ind_sort), 2, 2);    % take the second order differential to remove the background contributions
[~, snY] = estimate_baseline_sn(dY(:));  % estimate the noise level in dY 
dci = diff(ci(ind_sort), 2);
dci(dci>- thr_noise * sn) = 0;
ai = max(0, dY*dci'/(dci*dci'));  % use regression to estimate spatial component

% post-process ai by bwlabel
th = max(median(ai(:)), snY * 3 / sqrt(dci*dci'));  % minimum value of each pixel 
temp = full(ai>=th);
l = bwlabel(reshape(temp, nr, nc), 4);   % remove disconnected components
temp(l~=l(ind_ctr)) = false;
ai(~temp(:)) = 0;
if sum(ai(:)>0) < min_pixels %the ROI is too small
    ind_success=false;
    return;
end

% refine ci given ai 
ind_nonzero = (ai>0);
ai_mask = mean(ai(ind_nonzero))*ind_nonzero;
ci = (ai-ai_mask)'*ai\((ai-ai_mask)'*Y);
[b, sn] = estimate_baseline_nosie(ci); 
ci = ci - b; 
ind_neg = (ci<-4*sn); 
ci(ind_neg) = rand(sum(ind_neg), 1)*sn; 

% normalize the result
ci = ci / sn;
ai = ai * sn;
% return results
if norm(ai)==0
    ind_success= false;
else
    ind_success=true;
end