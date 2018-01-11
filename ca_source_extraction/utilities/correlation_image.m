function Cn = correlation_image(Y,sz,d1,d2,flag_norm, K)

% construct correlation image based on neighboing pixels
% Y: raw data
% sz: define the relative location of neighbours. it can be scalar, 2
%       element vector or a binary matrix
%       scalar: number of nearest neighbours, either 4 or 8, default 4
%       2 element vector: [rmin, rmax], the range of neighbours. the
%           distance between the neighbour and the pixel is d, dmin <= r <
%           dmax.
%       matrix: a squared matrix (2k+1)*(2k+1)indicating the location of neighbours
% d1,d2: spatial dimensions
% flag_norm: indicate whether Y has been normalized and centered ( 1 is
%   yes, 0 is no)
% K:  scalar, the rank of the random matrix for projection

% Author: Eftychios A. Pnevmatikakis, Simons Foundation, 2015
% with modifications from Pengcheng Zhou, Carnegie Mellon University, 2015.
% It uses convolution and random projection for speeding up the
% computation.

%% preprocess the raw data
if ~exist('flag_norm', 'var') || isempty(flag_norm)
    flag_norm = false;
end

if ~exist('sz', 'var') || isempty(sz)
    sz = [0,1,0; 1,0,1; 0,1,0];
end

% center data 
Y = bsxfun(@minus, double(Y), mean(Y, ndims(Y))); 
if ~ismatrix(Y)
    [d1, d2, T] = size(Y);
else
    T = size(Y, 2);
end

if exist('K', 'var') && (~isempty(K))
    Y = double(reshape(Y, [], T))*randn(T, K); 
    % centering
    mY = mean(Y,2);
    Y = bsxfun(@minus, Y, mY);        % normalizing
    flag_norm = false; 
end
if ~flag_norm
    sY = sqrt(mean(Y.*Y, ndims(Y)));
    sY(sY==0) = 1; % avoid nan values
    Y = bsxfun(@times, Y, 1./sY);
end
if ismatrix(Y)
    Y = reshape(Y, d1, d2, []);
end

%% construct a matrix indicating location of the matrix
if  isscalar(sz)
    if sz == 8      % 8 nearest neighbours
        sz = [1,1,1; 1,0,1; 1,1,1];
    elseif sz==4
        sz = [0,1,0; 1,0,1; 0,1,0];
    end
elseif length(sz(:)) == 2
    % the specified neighbours has a distance within the domain [dmin,
    % dmax)
    sz = ceil(sz);
    dmin = min(sz); dmax = max(sz);
    rsub = (-dmax+1):(dmax-1);      % row subscript
    csub = rsub;      % column subscript
    [cind, rind] = meshgrid(csub, rsub);
    R = sqrt(cind.^2+rind.^2);
    sz = (R>=dmin) .* (R<dmax);
end

%% compute the correlation
Yconv = imfilter(Y, sz);        % sum over the neighbouring pixels
MASK = imfilter(ones(d1,d2), sz);   % count the number of neighbouring pixels
Cn = mean(Yconv.*Y, 3)./MASK;   % compute correlation and normalize
