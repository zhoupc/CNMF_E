function [Yest, results, sn] = update_background(Y, varargin)
%% approximate the background with locally-linear embedding
%% inputs:
%   Y: d1*d2*T 3D matrix, video data
%   rr: scalar, average neuron size
%   ssub: spatial downsampling factor
%   ACTIVE_PX:  indicators of pixels to be approximated
%   sn:  noise level for each pixel
%   thresh: threshold for selecting frames with resting states

%% outputs:
%   Yest: d1*d2*T 3D matrix, reconstructed video data
%   results: struct variable {weights, ssub}
%       weights: d1*d2 cell, each element is a 2*J matrix. Row 1 has the indice of the
%       ring neighbors and row 2 has the corresponding weights.
%       ssub:    scalar, spatial downsampling factor

%% Author: Pengcheng Zhou, Carnegie Mellon University,2016

%% parse input arguments
[d1, d2, T] = size(Y);

% default values
Ybg = [];   % current estimation of the fluctuating background
rr = 15;
ssub = 1;
thresh = inf;
estimate_b0 = false;
sn = [];
b0 = mean(Y, 3);    % constant baseline for each pixel
k = 1;
while k<=length(varargin)
    switch lower(varargin{k})
        case 'ybg'  % previous estimation of the fluctuating background
            Ybg = varargin{k+1};
        case 'ssub'
            ssub = varargin{k+1};
        case  'rr'
            rr = varargin{k+1};
        case 'active_px'
            ACTIVE_PX = varargin{k+1};
        case 'sn'
            sn = varargin{k+1};
        case 'thresh'
            thresh = varargin{k+1};
        case 'estimate_b0'
            if k<length(varargin) && islogical(varargin{k+1})
                estimate_b0 = varargin{k+1};
            else
                estimate_b0 = true;
                k = k - 1;
            end
        otherwise
            k = k+1;
    end
    k = k+2;
end

%% detect outliers and replace it with the previous estimation
if thresh<inf    % do thresholding
    % estimate noise
    if isempty(sn)  % first iteration
        sn = get_noise_fft(reshape(Y, d1*d2,T));
    end
    sn = reshape(sn, d1, d2);
    
    %  current estimation
    if isempty(Ybg)  % use the mean of neighbors
        % construct the pattern of selecting neighbors
        rsub = (-rr):(rr);      % row subscript
        csub = rsub;      % column subscript
        [cind, rind] = meshgrid(csub, rsub);
        R = sqrt(cind.^2+rind.^2);
        neigh_kernel = (R>=rr) .* (R<rr+1);  % kernel representing the selected neighbors
        
        b0 = mean(Y, 3);
        Y = bsxfun(@minus, Y, b0);
        Ybf = bsxfun(@times, imfilter(Y, neigh_kernel), ...
            1./imfilter(ones(d1, d2), neigh_kernel));
        % detect events
        ind_events = (bsxfun(@times, Y-Ybf, 1./sn) > thresh);
        % replace the outliers with the current estimation
        Y(ind_events) = Ybf(ind_events); % remove potential calcium transients
        clear Ybf;
        
        temp = mean(Y, 3);
        Y = bsxfun(@minus, Y, temp);
        b0 = temp + b0;
    else
        % use previous estimation
        % detect events
        Ybg = reshape(Ybg, size(Y));
        ind_events = (bsxfun(@times, Y-Ybg, 1./sn) > thresh);
        Y(ind_events) = Ybg(ind_events);  % remove potential calcium transients
        clear Ybg; 
        b0 = mean(Y, 3);
        Y = bsxfun(@minus, Y, b0);
    end
else
    Y = bsxfun(@minus, Y, b0);
    ind_events = false(d1, d2, T);
end

%% downsampling
%downsample the data
Y0 = Y;
if ssub>1
    Y = imresize(Y, 1./ssub);
    [d1s, d2s, ~] = size(Y);
    rr = round(rr/ssub)+1;
    ind_events = (imresize(reshape(double(ind_events), d1, d2, []), [d1s, d2s])>0);
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

%% determine neibours of each pixel

% construct the pattern of selecting neighbors
rsub = (-rr):(rr);      % row subscript
csub = rsub;      % column subscript
[cind, rind] = meshgrid(csub, rsub);
R = sqrt(cind.^2+rind.^2);
neigh_kernel = (R>=rr) .* (R<rr+1);  % kernel representing the selected neighbors

[r_shift, c_shift] = find(neigh_kernel);
r_shift = r_shift - rr -1;
c_shift = c_shift - rr - 1;

[csub, rsub] = meshgrid(1:d2s, 1:d1s);
csub = reshape(csub, [], 1);
rsub = reshape(rsub, [], 1);
csub = bsxfun(@plus, csub, c_shift');
rsub = bsxfun(@plus, rsub, r_shift');
% remove neighbors that are out of boundary
ind = or(or(csub<1, csub>d2s), or(rsub<1, rsub>d1s));
csub(ind) = nan;
rsub(ind) = nan;

%% run approximation
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:SingularMatrix');

Y = reshape(Y, d1s*d2s, []);
ind_events = reshape(ind_events, d1s*d2s, []);
Yest = Y;
weights = cell(d1s, d2s);

for m=1:length(ind_px)
    % select pixel and its neighbors
    px = ind_px(m);
    ind_nhood = sub2ind([d1s,d2s], rsub(px, :), csub(px, :));
    ind_nhood(isnan(ind_nhood)) = [];
    
    % choose frames without events for estimating regression weights
    tmp_ind = ~ind_events(px, 2:end);
    X = Y(ind_nhood, tmp_ind);
    y = Y(px, tmp_ind);
    tmpXX = X*X';
    w = (tmpXX+eye(size(tmpXX))*sum(diag(tmpXX))*(1e-100)) \ (X*y');
    Yest(px, :) = w'*Y(ind_nhood, :);
    weights{px} = [ind_nhood; w'];
end
results.weights = weights;
results.ssub = ssub;

Yest = reshape(Yest, d1s, d2s, []);
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:SingularMatrix');

%% up-sampling
if ssub>1
    Yest = imresize(Yest, [d1, d2]);
    Y = Y0;
    clear Y0;
end

%% estimate the baseline and noise
if estimate_b0
    
    Y = reshape(Y, size(Yest));
    Yres = reshape(Y - Yest, d1*d2, T);  % residual
    clear Y; 
    b0_ds = zeros(d1s, d2s);
    sn_ds = zeros(d1s, d2s);
    indr = round(linspace(1,d1, d1s));
    indc = round(linspace(1,d2,d2s));
    ind = bsxfun(@plus, indr', (indc-1)*d1);
    parfor m=1:(d1s*d2s)
        % sub-sampling the data and estimate the noise
        [b0_ds(m), sn_ds(m)] = estimate_baseline_noise(Yres(ind(m), :));
    end
    b0 = b0 - imresize(b0_ds, [d1, d2]);
    sn = imresize(sn_ds, [d1, d2]);
end

%% return results
Yest = bsxfun(@plus, Yest, b0);














