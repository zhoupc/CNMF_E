function [c, s, b, g, active_set] = foopsi_oasisAR2(y, g, lam, smin, optimize_b,...
    optimize_g, decimate, maxIter)
%% Infer the most likely discretized spike train underlying an AR(2) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g_1 c_{t-1}- g_2 c_{t-2} >= 0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g:  1*2 vector, Parameter of the AR(2) process that models the fluorescence ...
%impulse response.
%   lam:  scalar, sparsity penalty parameter
%   smin: scalar, minimum spike size 
%   optimize_b: bool, optimize baseline if True
%   optimize_g: integer, number of large, isolated events to consider for
%       optimizing g
%   decimate: int, decimation factor for estimating hyper-parameters faster
%       on decimated data
%   maxIter:  int, maximum number of iterations
%   active_set: npool x 4 matrix, warm stared active sets

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   b: scalar, fluorescence baseline
%   g: scalar, parameter of the AR(1) process
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging


%% input arguments
y = reshape(y, [], 1);
T = length(y);

if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y, 2);
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end

if ~exist('smin', 'var') || isempty(smin);   
    smin = 0; 
end
if ~exist('optimize_b', 'var') || isempty(optimize_b)
    optimize_b = false;
end
if ~exist('optimize_g', 'var') || isempty(optimize_g)
    optimize_g = 0;
end
if ~exist('decimate', 'var') || isempty(decimate)
    decimate = 1;
else
    decimate = max(1, round(decimate));
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 10;
end

% change parameters due to downsampling
if decimate>1
    decimate = 1;  %#ok<NASGU>
    disp('to be done');
    %     fluo = y;
    %     y = resample(y, 1, decimate);
    %     g = g^decimate;
    %     thresh = thresh / decimate / decimate;
    %     T = length(y);
end

%% optimize parameters
if ~optimize_b   %% don't optimize the baseline b
    %% initialization
    b = 0;
    [solution, spks, active_set] = oasisAR2(y, g, lam, smin);
    
    %% iteratively update parameters g
    if optimize_g;
        g_old = g;
        g = update_g(y-b, g, spks);
        [solution, spks,active_set] = oasisAR2(y, g, 0, smin);
    end
else
    disp('to be done'); 
%     %% initialization
%     b = quantile(y, 0.15);
%     [solution, spks, active_set] = oasisAR2(y-b, g, lam, smin);
%     
    %% iteratively update parameters g and b
%     for m=1:maxIter
%         b = mean(y-solution);
%         if and(optimize_g, ~g_converged);
%             g0 = g;
%             g = update_g(y-b, g, spks);
%             smin = choose_smin(g, sn, 0.9999);
%             [solution, spks,active_set] = oasisAR2(y, g, 0, smin);
%             
%             if norm(g-g0,2)/norm(g0) < 1e-3 % g is converged
%                 g_converged = true;
%             end
%         else
%             break;
%         end
%     end
end
c = solution;
s = spks;
end


%update the AR coefficient: g
function g = update_g(y, g, spks)
%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g: 2x1 vector, AR2 parameters
%   active_set: npools*4 matrix, previous active sets.
% smin: scalr, minimize size of nonzero spikes

%% outputs
%   c: T*1 vector
%   s: T*1 vector, spike train
%   active_set: npool x 4 matrix, active sets
%   g: scalar

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016

%% initialization
s_th = quantile(spks(spks>1e-3), 0.25);
tsp = find(spks>=s_th);
tsp = reshape(tsp, 1, []);
time_p = find(conv2(double(spks<=s_th), ones(30,1), 'same')>0);
time_p = reshape(time_p,[],1);
y = reshape(y,[],1);    % fluorescence data
yp = y(time_p);
T = length(y);
tau_dr = ar2exp(g);
tau_d = tau_dr(1);
tau_r = tau_dr(2);
warning('off', 'MATLAB:singularMatrix');

%% find the optimal g and get the warm started active_set
tau_d0 = tau_d;
tau_r0 = tau_r;
bnd_d = tau_d0 * [1/4, 4];
bnd_r = tau_r0 * [1/4, 4];
for m=1:10
    tau_r = fminbnd(@rss_taur, bnd_r(1), bnd_r(2));
    tau_d = fminbnd(@rss_taud, bnd_d(1), bnd_d(2));
    if and(abs(tau_d-tau_d0)/tau_d0 < 1e-4, abs(tau_r-tau_r0)/tau_r0 < 1e-4)
        break;
    else
        tau_d0 = tau_d;
        tau_r0 = tau_r;
    end
end

%% compute the optimal solution
g = exp2ar([tau_d, tau_r]);

%% nested functions

    function rss = rss_taur(tau_r)
        ht = (exp(-(1:T)/tau_d) - exp(-(1:T)/tau_r))/(tau_d-tau_r);
        ht(T) = 0;
        ind = bsxfun(@minus, time_p, tsp);
        ind(ind<=0) = T;
        V = ht(ind);
        
        % find the best value of st
        s = (V'*V)\(V'*yp);
        res = yp - V*s;
        rss = res' * res;
    end

    function rss = rss_taud(tau_d)
        ht = (exp(-(1:T)/tau_d) - exp(-(1:T)/tau_r))/(tau_d-tau_r);
        ht(T) = 0;
        ind = bsxfun(@minus, time_p, tsp);
        ind(ind<=0) = T;
        V = ht(ind);
        
        % find the best value of st
        s = (V'*V)\(V'*yp);
        res = yp - V*s;
        rss = res' * res;
    end
end