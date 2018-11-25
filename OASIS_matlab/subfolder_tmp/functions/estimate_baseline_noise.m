function [b, sn] = estimate_baseline_noise(y, bmin)
%% estimate the baseline and noise level for single calcium trace 
%% inputs: 
%   y:  T*1 vector, calcium trace 
%   thr: scalar, threshold for choosing points for fitting a gaussian
%   distribution. 
%   bmin: scalar, minimum value of the baseline, default(-inf)
%% outputs: 
%   b: scalar, baseline 
%   sn: scalar, sigma of the noise 
%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016 

%% input arguments 
if ~exist('bmin', 'var') || isempty(bmin) 
    bmin = -inf; 
end

%% create the histogram for fitting 
temp = quantile(y, 0:0.1:1); 
dbin = max(min(diff(temp))/3, (max(temp)-min(temp))/1000); 
bins = temp(1):dbin:temp(end); 
nums = hist(y, bins); 

%% fit a gaussian distribution: nums = A * exp(-(bins-b)/(2*sig^2))
if isempty(bins)
    b = mean(y);
    sn = 0;
    return; 
end
[b, sn] = fit_gauss1(bins, nums, 0.3, 3);

if b<bmin
    b = bmin;
    sn = fit_gauss1(bins-bmin, nums, 0.3, 3,false );
end
