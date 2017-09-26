function [W, b0] = fit_ring_model(Y, A, C, W_old, thresh_outlier, sn, ind_patch)
%% fit a ring model Bf = W*Bf s.t. Bf = Y-AC-b0*1^T to reconstruct the background 
%% inputs: 
%   Y:  d*T matrix 
%   A:  d*K matrix 
%   C:  d*T matrix 
%   W:  d*d matrix 
%   b0: d*1 vector 
%   thresh_outlier 
%   sn: d*1 noise level 

%% compute the fluctuating background 
Ymean = mean(Y,2); 
Cmean = mean(C, 2); 
Y = bsxfun(@minus, double(Y), Ymean); 
C = bsxfun(@minus, C, Cmean);
Bf = Y - A*C; 
b0 = Ymean(ind_patch) - A(ind_patch,:)*Cmean;

%% compute the previous estimatin and take care of the outliers. 
Bf_old = W_old*Bf; 
tmp_Bf = Bf(ind_patch, :); 
ind_outlier = bsxfun(@gt, tmp_Bf, bsxfun(@plus, Bf_old, thresh_outlier*reshape(sn, [], 1))); 
tmp_Bf(ind_outlier) = Bf_old(ind_outlier);  
Bf(ind_patch, :) = tmp_Bf; 
clear tmp_Bf Bf_old ind_outlier; 

%% fit a regression model to get W 
T = size(Y, 2); 
ind_pixels = find(ind_patch); 
d = length(ind_pixels); 
vec_ones = ones(1, T); 
W = W_old; 

for m=1:d
    idx = ind_pixels(m); 
    ind_ring = (W_old(m,:)~=0); 
    y = Bf(idx, :); 
    X = [Bf(ind_ring,:); vec_ones]; 
    
    tmpXX = X*X'; 
    tmpXy = X*y'; 

    w = (tmpXX+eye(size(tmpXX))*sum(diag(tmpXX))*(1e-5)) \ tmpXy;
    W(m, ind_ring) = w(1:(end-1))+(1e-100); % add a small value to keep those nonzero elements still nonzero
end 
%% 

