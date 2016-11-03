function [mu, sig, A] = fit_gauss1(x, y, thr, maxIter, mu_fix)
%% fit a gaussin curve given the data (x, y): y = A*exp(-(x-mu)^2/2/sig^2)

%% inputs:
%   x: T*1 vector
%   y: T*1 vector
%   thr: scalar, threshold for selecting points for fitting
%   maxIter: scalar, maximum iteration
%   mu_fix: boolean, fix the center to be 0 or not
%% outputs:
%   mu: scalar, center of the peak
%   sig: scalar, gaussian width
%   A:  scalar, amplitude of the gaussian curve

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
%% Reference:
%   A Simple Algorithm for Fitting a Gaussian Function , Hongwei Guo, 2011

%% input arguments
x = reshape(x, [], 1);
y = reshape(y, [], 1);
T = length(y);

if ~exist('thr', 'var') || isempty(thr)
    thr = 0.1;
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 5;
end
if ~exist('mu_fix', 'var') || isempty(mu_fix)
    mu_fix = false;
end
ind = (y>max(y)*thr);
x = x(ind);
y = y(ind);

x2 = x.^2;
x3 = x.^3;
x4 = x.^4;
y2 = y.^2;
logy = log(y);
y2logy = y2.*logy;
vec1 = ones(1, length(y));

%% fit the curve
if mu_fix % fix the mu to be 0 
    for miter=1:maxIter
        M = [vec1*y2, x2'*y2; ...
            x2'*y2, x4'*y2];
        b = [vec1*y2logy; x2'*y2logy];
        p = M\b;
        
        logy = p(1)*vec1' + p(2)*x2;
        y = exp(logy);
        y2 = y.^2;
        y2logy = y2.*logy;
    end
    
    mu= 0;
    sig = sqrt(-0.5/p(2));
    A = exp(p(1));
else
    for miter=1:maxIter
        M = [vec1*y2, x'*y2, x2'*y2; ...
            x'*y2, x2'*y2, x3'*y2; ...
            x2'*y2, x3'*y2, x4'*y2];
        b = [vec1*y2logy; x'*y2logy; x2'*y2logy];
        p = M\b;
        
        logy = p(1)*vec1' + p(2)*x + p(3)*x2;
        y = exp(logy);
        y2 = y.^2;
        y2logy = y2.*logy;
    end
    mu= -p(2)/2/p(3);
    sig = abs(sqrt(-0.5/p(3)));
    A = exp(p(1)-0.25*p(2)^2/p(3));
end
