function [C, C_raw, results_deconv] = HALS_temporal(Y, A, C, maxIter, deconv_options)
%% run HALS by fixating all spatial components 
% input: 
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components 
%   C:  K*T, temporal components 
%   maxIter:    maximum iteration number 
%   deconv_options:  struct variable, options for running deconvolution 
% output: 
%   C: K*T, updated temporal components 
%   C_raw: K*T, the temporal components before thresholded or being
%   denoised. 
%   results_deconv: struct variable with following fields: 
%       S: K*T, spiking activity 
%       sn: K*1 vector, noise level for each neurons 
%       kernel_pars: K*1 cell, parameters. 

% Author: Pengcheng Zhou, Columbia University, 2017 
% zhoupc1988@gmail.com 

%% options 
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter=1;
end
if ~exist('deconv_options', 'var') || isempty(deconv_options)
    % no deconvolution 
    deconv_flag = false;
else
    deconv_flag = true; 
end 
%% initialization 
K = size(A, 2);     % number of components 
T = size(Y, 2);     % number of frames 
if ~exist('C', 'var')||isempty(C)
    C = zeros(K, T);
end
C_raw = zeros(K,T);  % raw temporal components before being denoised and deconvolved.

A = full(A);       % spatial shapes of all neurons 
U = A'*Y;          
V = A'*A;
aa = diag(V);       % squares of l2 norm all all components
ind_update = find(aa>0);  % avoid updating neurons with empty spatial shapes. This only happens in very weird cases. 
n_update = length(ind_update);   % number of neurons to be updated. 
if deconv_flag
    S = zeros(K, T);
    sn = zeros(1, K); 
    kernel_pars = cell(1,K); 
end
%% updating 
for miter=1:maxIter
    for m=1:n_update
        k = ind_update(m); 
        ck_raw = C(k, :) + (U(k, :)-V(k, :)*C)/aa(k);     % update ck
                
        if ~deconv_flag
            % no deconvolution, just thresholding
            C_raw(k, :) = ck_raw;
            C(k, :) = max(0, ck_raw);
            if miter==maxIter
                C_raw(k, :) = ck_raw;
            end
        else
            %remove baseline and estimate noise level. estiamte the noise
            %using two methods: psd and histogram. choose the one with
            %smaller value. 
            try
                [b_hist, sn_hist] = estimate_baseline_noise(ck_raw);
            catch
                sn_hist = inf;
            end
            b = mean(ck_raw(ck_raw<median(ck_raw)));
            sn_psd = GetSn(ck_raw);
            if sn_psd<sn_hist
                tmp_sn = sn_psd;
            else
                tmp_sn = sn_hist;
                b = b_hist;
            end
            
            % subtract the baseline 
            ck_raw = ck_raw -b;
            sn(k) = tmp_sn;
            
            % deconvolution
            [ck, sk, tmp_options]= deconvolveCa(ck_raw, deconv_options, 'maxIter', 2, 'sn', tmp_sn, 'pars', kernel_pars{k});
            kernel_pars{k} = reshape(tmp_options.pars, 1, []);
            ck_raw = ck_raw - tmp_options.b;
            if sum(abs(ck))==0 % avoid the case where neuron's temporal activity is 0
                ck = ck_raw; 
            end 
            C(k, :) = ck;         
            % save the spike count in the last iteration
            if miter==maxIter
                S(k, :) = sk;
                C_raw(k, :) = ck_raw;
            end

        end
    end
end

if deconv_flag
    results_deconv.sn = sn; 
    results_deconv.kernel_pars = kernel_pars; 
else
    results_deconv = []; 
end 
