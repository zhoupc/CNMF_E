function [W, b0] = fit_ring_model(Y, A, C, W_old, thresh_outlier, sn, ind_patch, with_projection)
%% fit a ring model Bf = W*Bf s.t. Bf = Y-AC-b0*1^T to reconstruct the background 
%% inputs: 
%   Y:  d*T matrix 
%   A:  d*K matrix 
%   C:  d*T matrix 
%   W:  d*d matrix 
%   b0: d*1 vector 
%   thresh_outlier 
%   sn: d*1 noise level 
%   with_projection: boolean 

%% find pixels to be updated 
[d, T] = size(Y); 
if ~exist('A', 'var') || isempty(A)
    A = ones(d,1); 
    C = zeros(1, T); 
elseif issparse(A)
    A = full(A); 
end
if issparse(A)
    A = full(A); 
end
% check whether this is the first run 
if length(unique(W_old(1,:)))==2
    ind_active = true(size(W_old(:,1))); 
else
    ind_active = (abs(W_old)*sum(A,2)>0);
end

if ~exist('sn', 'var') || isempty(sn)
    sn = GetSn(double(Y)); 
end
if ~exist('ind_patch', 'var') || isempty(ind_patch)
    ind_patch = true(size(Y,1), 1); 
end
if ~exist('with_projection', 'var') || isempty(with_projection)
    with_projection = true; 
end

%% compute the fluctuating background 
Ymean = mean(Y,2); 
Cmean = mean(C, 2); 
b0 = Ymean(ind_patch) - A(ind_patch,:)*Cmean;
Y = bsxfun(@minus, double(Y), Ymean); 
C = bsxfun(@minus, C, Cmean);
Bf = Y - A*C; 

%% compute the previous estimation and take care of the outliers. 
if ~isnan(thresh_outlier)
    Bf_old = W_old*Bf;
    tmp_Bf = Bf(ind_patch, :);
    ind_outlier = bsxfun(@gt, tmp_Bf, bsxfun(@plus, Bf_old, thresh_outlier*reshape(sn, [], 1)));
    tmp_Bf(ind_outlier) = Bf_old(ind_outlier);
    Bf(ind_patch, :) = tmp_Bf;
end

%% we don't need all frames for weights estimation.
T = size(Y, 2);
pmax = max(sum(W_old>0, 2));
nmax = pmax*100;
if (~isnan(thresh_outlier)) && (nmax<T)
    temp = sum(ind_outlier);
    ind_frames = (temp<=quantile(temp, nmax/T));
    nmax = nnz(ind_frames);
    Bf = Bf(:, ind_frames);
    vec_ones = ones(1, nmax);
else
    vec_ones = ones(1, T);
end
ind_pixels = find(ind_patch);
d = length(ind_pixels);
W = W_old;
T = size(Bf, 2); 
clear tmp_Bf Bf_old;

%% with prejection 
if with_projection
%     % subsampling data 
%     ind = linspace(1,T , min(nmax, T/2)); 
%     Bf = Bf(:, ind_k); 
    %     V = randn(size(Bf, 2), nk);
    %     Bf_proj = Bf*V;
    nk = min(round(T/1), nmax); 
    k = floor(T/nk);
    if k~=1
        Bf = Bf(:, 1:k:end);
%         Bf = imresize(Bf, [size(Bf, 1), nk], 'nearest'); 
%         Bf = squeeze( mean(reshape(Bf(:, 1:(k*nk)), [], k, nk), 2)); %imresize(Bf, [size(Bf, 1), nk]);
    end
    vec_ones = ones(1, size(Bf, 2));
    for m=1:d
        if ~ind_active(m)
            continue;
        end
        idx = ind_pixels(m);
        % choose neighbors with low neural activity
        
        ind_ring = (W_old(m,:)~=0);
        y = Bf(idx, :);
        X = [Bf(ind_ring,:); vec_ones];
        
        tmpXX = X*X';
        tmpXy = X*y';
        
        w = (tmpXX+eye(size(tmpXX))*sum(diag(tmpXX))*(1e-5)) \ tmpXy;
        W(m, ind_ring) = w(1:(end-1))+(1e-100); % add a small value to keep those nonzero elements still nonzero
    end
else
    %% fit a regression model to get W
    
    for m=1:d
        if ~ind_active(m)
            continue;
        end
        idx = ind_pixels(m);
        ind_ring = (W_old(m,:)~=0);
        y = Bf(idx, :);
        X = [Bf(ind_ring,:); vec_ones];
        
        tmpXX = X*X';
        tmpXy = X*y';
        
        w = (tmpXX+eye(size(tmpXX))*sum(diag(tmpXX))*(1e-5)) \ tmpXy;
        W(m, ind_ring) = w(1:(end-1))+(1e-100); % add a small value to keep those nonzero elements still nonzero
    end
end
%% 

