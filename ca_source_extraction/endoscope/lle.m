function Yest = lle(Y, ssub, rr, ACTIVE_PX, method)
%% approximate the background with locally-linear embedding
%% inputs:
%   Y: d1*d2*T 3D matrix, video data
%   rr: scalar, average neuron size
%   ssub: spatial downsampling factor
%   ACTIVE_PX:  indicators of pixels to be approximated
%% outputs:
%   Yest: d1*d2*T 3D matrix, reconstructed video data
%% Author: Pengcheng Zhou, Carnegie Mellon University,2016

%% input arguments
[d1, d2, T] = size(Y);

% center the fluorescence intensity by its mean 
Ymean = mean(Y, 3); 
Y = Y - bsxfun(@minus, Ymean, ones(1, 1, T)); 

% average neuron size
if ~exist('rr', 'var')|| isempty(rr)
    rr = 15;
end
% spatial downsampling
if ~exist('ssub', 'var') || isempty(ssub)
    ssub = 1;
end

%downsample the data
if ssub>1
    Y = imresize(Y, 1./ssub);
    [d1s, d2s, ~] = size(Y);
    rr = round(rr/ssub)+1;
else
    d1s = d1;
    d2s = d2;
end

% pixels to be approximated
if exist('ACTIVE_PX', 'var') && ~isempty(ACTIVE_PX)
    ACTIVE_PX = reshape(double(ACTIVE_PX), d1, d2);
    ACTIVE_PX = (imresize(ACTIVE_PX, 1/ssub)>0);
    ind_px = find(ACTIVE_PX(:));  % pixels to be approxiamted
else
    ind_px = (1:(d1s*d2s));
end

% select method for approximation 
if ~exist('method', 'var') || isempty(method)
    method = 'regression'; 
end 

%% determine neibours of each pixel
c_shift = (-rr:rr);
r_shift = round(sqrt(rr^2-c_shift.^2));
c_shift = [c_shift, c_shift(2:end)];
r_shift = [-r_shift, r_shift(2:end)];

[csub, rsub] = meshgrid(1:d2s, 1:d1s);
csub = reshape(csub, [], 1);
rsub = reshape(rsub, [], 1);
csub = bsxfun(@plus, csub, c_shift);
rsub = bsxfun(@plus, rsub, r_shift);
% remove neighbors that are out of boundary
ind = or(or(csub<1, csub>d2s), or(rsub<1, rsub>d1s));
csub(ind) = nan;
rsub(ind) = nan;

%% run approximation
if strcmpi(method, 'mean')
    mean_kernel = zeros(length(c_shift));
    mean_kernel(sub2ind(size(mean_kernel), r_shift+rr+1, c_shift+rr+1)) = 1/length(c_shift);
    Yest = imfilter(Y, mean_kernel, 'replicate');
    Yest = reshape(Yest, d1s, d2s, []); 
else
    gamma = 0.001; % add regularization
    Y = reshape(Y, d1s*d2s, []);
    Yest = zeros(size(Y));
    for m=1:length(ind_px)
        px = ind_px(m);
        ind_nhood = sub2ind([d1s,d2s], rsub(px, :), csub(px, :));
        ind_nhood(isnan(ind_nhood)) = [];
        J = length(ind_nhood);
        
        if strcmpi(method, 'regression')
            G = bsxfun(@minus, Y(ind_nhood, 1:3:end), Y(px, 1:3:end));
            w = (G*G'+gamma*sum(G(:).^2)*eye(J))\ones(length(ind_nhood), 1);
            w = w/sum(w);
            Yest(px, :) = w'*Y(ind_nhood, :);
        else%if strcmpi(method, 'median')
            Yest(px, :) = median(Y(ind_nhood, :), 1);
        end
        
    end
    ind = 1:(d1s*d2s);
    ind(ind_px) = [];
    if ~isempty(ind)
        temp = imfilter(Y, ones(3, 3)/9, 'replicate');
        Yest(ind, :) = temp(ind, :); % without approximation
    end
    Yest = reshape(Yest, d1s, d2s, []);
end
%% return the result
if ssub>1 %up sampling
    Yest = imresize(Yest, [d1, d2]);
end
Yest = Yest + bsxfun(@times, Ymean-median(Yest, 3), ones(1, 1, T)); 