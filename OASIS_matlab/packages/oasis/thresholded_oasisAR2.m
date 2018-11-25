function [c, s, b, g, smin, active_set] = thresholded_oasisAR2(y, g, sn, smin, optimize_b,...
    optimize_g, decimate, maxIter, thresh_factor)
%% Infer the most likely discretized spike train underlying an AR(2) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g_1*c_{t-1} - g_2*c_{t-2}>=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   sn:  scalar, standard deviation of the noise distribution
%   optimize_b: bool, optimize baseline if True
%   optimize_g: integer, number of large, isolated events to consider for
%       optimizing g
%   decimate: int, decimation factor for estimating hyper-parameters faster
%       on decimated data
%   maxIter:  int, maximum number of iterations
%   active_set: npool x 4 matrix, warm stared active sets
%  thresh_factor: scalar, set the maximum thresh as thresh_factor*sn^2*T

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   b: scalar, fluorescence baseline
%   g: scalar, parameter of the AR(1) process
%   smin: scalar, minimum nonzero spike count
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

if ~exist('optimize_b', 'var') || isempty(optimize_b)
    optimize_b = false;
end

if ~exist('sn', 'var') || isempty(sn)
    if optimize_b
        [b, sn] = estimate_baseline_noise(y);
    else
        [~, sn] = estimate_baseline_noise(y-smooth(y,100));
    end
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
if ~exist('thresh_factor', 'var') || isempty(thresh_factor)
    thresh_factor = 1.0;
end

%% start from smin that avoid counting gaussian noise as a spike
smin = choose_smin(g, sn, 0.99999999);
thresh = thresh_factor* sn * sn * T;

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
g_converged = false;

%% optimize parameters
tol = 1e-4;
if ~optimize_b   %% don't optimize the baseline b
    %% initialization
    b = 0;
    [solution, spks, active_set] = oasisAR2(y, g, [], smin);
    
    res = y - solution;
    RSS0 = res' * res;
    %% iteratively update parameters lambda & g
    for miter=1:maxIter
        if isempty(active_set)
            break;
        end
        % update g
        if miter==1 && and(optimize_g, ~g_converged);
            g0 = g;
            g = update_g(y-b, g, spks);
            smin = choose_smin(g, sn, 0.9999);
            [solution, spks,active_set] = oasisAR2(y, g, 0, smin);
            if norm(g-g0,2)/norm(g0) < 1e-3 % g is converged
                g_converged = true;
            end
        end
        
        res = y - solution;
        RSS = res' * res;
        if abs(RSS-RSS0)<tol  % no more changes
            break;
        end
        len_active_set = size(active_set,1);
        
        if or(RSS>thresh, sum(solution)<1e-9)  % constrained form has been found, stop
            break;
        else
            RSS0 = RSS;
            % update smin
            [smin, solution, spks, active_set] = update_smin(y, g, smin,...
                solution, spks, active_set, sqrt(thresh), max(spks));
        end
    end
else
    %% initialization
    [b, ~] = estimate_baseline_noise(y); 
    [solution, spks, active_set] = oasisAR2(y-b, g, [], smin);
    
    res = y - solution -b;
    RSS0 = res' * res;
    %% iteratively update parameters lambda & g
    for miter=1:maxIter
        if isempty(active_set)
            break;
        end
        % update g
        if and(optimize_g, ~g_converged);
            g0 = g;
            g = update_g(y-b, g, spks);
            smin = choose_smin(g, sn, 0.9999);
            [solution, spks,active_set] = oasisAR2(y, g, 0, smin);
            
            if norm(g-g0,2)/norm(g0) < 1e-3 % g is converged
                g_converged = true;
            end
        end
        
        res = y - solution -b;
        RSS = res' * res;
        if abs(RSS-RSS0)<tol  % no more changes
            break;
        end
        len_active_set = size(active_set,1);
        
        if or(RSS>thresh, sum(solution)<1e-9)  % constrained form has been found, stop
            break;
        else
            RSS0 = RSS;
            % update smin
            [smin, solution, spks, active_set] = update_smin(y-b, g, smin,...
                solution, spks, active_set, sqrt(thresh), max(spks));
            b = mean(y-solution);
        end
    end
end
c = solution;
s = spks;
g = g(1:2);

%% nested functions
    function [smin, solution, spks, active_set] = update_smin(y, g, smin, solution, ...
            spks, active_set, thresh, s_max)
        %%estimate smin to match the thresholded RSS
        len_active_set = size(active_set, 1);
        sv = linspace(smin, s_max, min(9, len_active_set));
        ind_start = 1;
        ind_end = length(sv);
        
        while (ind_end-ind_start)>1
            ind = floor((ind_start+ind_end)/2);
            tmp_smin = sv(ind);
            [tmp_solution, tmp_spks, tmp_active_set] = oasisAR2(y, g, [], ...
                tmp_smin, [], [], active_set);
            sqRSS = norm(y-tmp_solution,2);
            if sqRSS<thresh % increase smin
                solution = tmp_solution;
                spks = tmp_spks;
                active_set = tmp_active_set;
                smin = tmp_smin;
                ind_start = ind;
            elseif sqRSS>thresh % decrease smin
                ind_end = ind;
            else
                break;
            end
        end
    end


end
%
% function [c, active_set, g, s] = update_g(y, g, spks, smin, c)
% %% update the AR coefficient: g
%
%
% %% residual
% yres = y - c;
% T = length(y);
% %% convolution kernel
% ht = filter(1,[1,-g],[1,zeros(1,500)]);
% ht(ht<0.01) = [];
% w = length(ht);
% ht = [ zeros(1, 2), ht];
% %% find all spikes
% tsp = find(spks>0);
% tsp(tsp<3) = [];
% tsp(tsp==T) = [];
% spv = spks(tsp);
%
%
% %% compute the mean waveform
% mean_trace = zeros(1, 2+w);
% for m=1:length(tsp)
%     ti = tsp(m);
%     if ti<= T-w
%         mean_trace = mean_trace + ht*spv(m) + yres((ti-2):(ti+w-1))' ;
%     else
%         ind = (ti-2):T;
%         mean_trace(1:(T-ti+3)) = mean_trace(1:(T-ti+3)) + ht(1:(T-ti+3))*spv(m) + yres(ind)'/spv(m);
%     end
% end
%
% g = estimate_time_constant(mean_trace);
% [c, s, active_set] = oasisAR2(y, g, [], smin);
% end

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
    try
        tau_r = fminbnd(@rss_taur, bnd_r(1), bnd_r(2));
        tau_d = fminbnd(@rss_taud, bnd_d(1), bnd_d(2));
    catch
        break;
    end
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