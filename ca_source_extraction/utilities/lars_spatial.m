function A = lars_spatial(Y, A, C, active_pixel, sn)
%% run HALS by fixating all spatial components 
% input: 
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components 
%   C:  K*T, temporal components 
%   active_pixel, mask for pixels to be updated 
%   sn: noise level for each pixel 

% output: 
%   A: d*K, updated spatial components 

% Written by:
% Pengcheng Zhou, Columbia University, 2017

%%
[d, T] = size(Y); 
K = size(C, 1); 
if isempty(A)
    A = zeros(d, K); 
end
if isempty(C)
    A = zeros(d, 0); 
    return; 
end
if nargin<5
    % noise level 
    sn = GetSn(double(Y));
end   
sn = reshape(sn, [], 1); 

if nargin<4
    active_pixel=true(d, T);
elseif isempty(active_pixel)
    active_pixel = true(d, T); 
else
    active_pixel = logical(active_pixel); 
end     %determine nonzero pixels
 
K = size(C, 1);     % number of components 
[d, T] = size(Y); 

%% run nnls to get solution that preserves the sparsity of A 
Ymean = mean(Y,2); 
Y = bsxfun(@minus, Y, Ymean); 
C = bsxfun(@minus, C, mean(C,2)); 
CC = C*C';
YC = C*Y';
ind_fit = find(sum(active_pixel,2)>1e-9);
A = zeros(d, K); 
thresh = sn.^2*T - sum(Y.^2, 2); 

for m=1:length(ind_fit)
    ind = active_pixel(ind_fit(m), :);
    A(ind_fit(m), ind) = nnls(CC(ind, ind), YC(ind, ind_fit(m)), [], [], [], thresh(m));
end


function s = nnls(A, b, s, tol, maxIter, thresh)
%% fast algorithm for solving nonnegativity constrained least squared
% problem minize sum(s) s.b. norm(y-K*s, 2)<=sn^2*T.

%% inputs:
%   A: n x p matrix, K'*K
%   b: n x 1 vector, K'*y
%   s: p x 1 vector, warm started s
%   tol: scalar, smallest nonzero values
%   maxIter: scalar, maximum nonzero values

%% outputs:
%   s: p x 1 vector, solution

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging
% Bro R & Jong S, Journal of Chemometrics 1997, A FAST NON-NEGATIVITY-CONSTRAINED LEAST SQUARES ALGORITHM

%% input arguments
p = size(A,2);      % dimension of s
if ~exist('s', 'var') || isempty(s)
    s = zeros(p, 1);
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-9;
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = p;
end
if ~exist('thresh', 'var') || isempty(thresh)
    thresh = []; 
end
if sum(s>0)>maxIter
    s = zeros(p,1);
end
for miter=1:maxIter
    l = b - A*s;            % negative gradient
    Pset = (s>0);       % passive set
    
    if max(l) < tol         % no more passive set
        break;
    end
    
    if ~isempty(thresh)
        if s'*A*s-2*s'*b<=thresh
            break;
        end
    end
    [~, temp] = max(l);         % choose the one with the largest gradient
    Pset(temp) = true;         % add it to the passive set
    if sum(Pset)>maxIter
        break;
    end
    % correct nonnegativity violations
    while any(Pset)
        % run unconstrained least squares for variables in passive sets
        try
            mu = A(Pset, Pset) \ b(Pset);
        catch
            mu = (A(Pset, Pset) + tol*eye(sum(Pset))) \ b(Pset);
        end
        
        if all(mu>tol)
            break;
        end
        
        temp = s(Pset) ./ (s(Pset)-mu);
        temp(mu>tol) = [];
        a = min(temp);
        s(Pset) = s(Pset) + a*(mu-s(Pset));
        Pset(s<tol) = false;
    end
    
    s(Pset) = mu;
end















