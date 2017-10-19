function [Ydt, X, R] = detrend_data(Y,nk, method)
%% detrend fluorescence signals with B-spline basis 
%% Inputs: 
%   Y: d X T matrix, video data 
%   nk: scalar, number of knots
%   method: {'local_min', 'spline'} method for removing the trend 
% Outputs: 
%   Ydt: d X T matrix, detrended data 
%   X: T*M matrix, each column is one basis 
%   R: d*M matrix, coefficients of all basis for all pixels 

%% create basis 
[~, T] = size(Y); 

if ~exist('nk', 'var')
    nk = 5; 
end 
if ~exist('method', 'var') || isempty(method) 
    method = 'spline'; 
end

if strcmpi(method, 'spline')
    X = bsplineM((1:T)', linspace(1, T, nk), 4);
    
    %% compute coefficients of all spline basis
    R = (Y*X)/(X'*X);
    
    %% compute detrended data
    Ydt = Y-R*X';
else
    k = ceil(T/nk);
    [d, T] = size(Y);
    Tnew = ceil(T/k)*k;
    if T~=Tnew
        Y(:, (T+1):Tnew) = repmat(Y(:, T), [1, Tnew-T]);
    end
    Y = reshape(Y, d, k, []);
    Ydt = reshape(bsxfun(@minus, Y, min(Y,[], 2)), d, []);
    Ydt = Ydt(:, 1:T);
    X = [];
    R = [];
end 
