function A = HALS_spatial(Y, A, C, active_pixel, maxIter)
%% run HALS by fixating all temporal components 
% input: 
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components 
%   C:  K*T, temporal components 
%   active_pixel, mask for pixels to be updated 

% output: 
%   A: d*K, updated spatial components 

% Author: Pengcheng Zhou, Carnegie Mellon University, adapted from Johannes
% Friedrich's NIPS paper "Fast Constrained Non-negative Matrix
% Factorization for Whole-Brain Calcium Imaging Data

%% options for HALS
if nargin<5;    maxIter = 1;    end   %maximum iteration number 
if nargin<4;    active_pixel=true(size(A));
elseif isempty(active_pixel)
    active_pixel = true(size(A)); 
else
    active_pixel = logical(active_pixel); 
end     %determine nonzero pixels 

%% initialization 
A(~active_pixel) = 0; 
K = size(A, 2);     % number of components 
Cmean = mean(C,2); 
Ymean = mean(Y,2); 
T = size(C,2); 
U = double(Y*C'-T*(Ymean*Cmean')); 
V = double(C*C'-T*(Cmean*Cmean')); 
cc = diag(V);   % squares of l2 norm all all components 

%% updating 
for miter=1:maxIter
    for k=1:K
        if cc(k)==0
            continue; 
        end
        tmp_ind = active_pixel(:, k); 
        ak = max(0, A(tmp_ind, k)+(U(tmp_ind, k)-A(tmp_ind,:)*V(:, k))/cc(k)); 
        A(tmp_ind, k) = ak; 
    end
end