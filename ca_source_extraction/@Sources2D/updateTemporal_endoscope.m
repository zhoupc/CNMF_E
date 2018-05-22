function [C_offset,ind_del] = updateTemporal_endoscope(obj, Y, allow_deletion)
%% run HALS by fixating all spatial components
% input:
%   Y:  d*T, fluorescence data
%   allow_deletion: boolean, allow deletion (default: true)
% output:
%   C_raw: K*T, temporal components without being deconvolved

% Author: Pengcheng Zhou, Carnegie Mellon University, adapted from Johannes

% options
global mode nam
if strcmp(mode,'initiation')
    maxIter = obj.options.maxIter;
else
    maxIter=1;
end
deconv_options_0 = obj.options.deconv_options;

if ~exist('allow_deletion', 'var')
    allow_deletion = false;
end
%% initialization
A = obj.A;
K = size(A, 2);     % number of components
C = obj.C;
C_raw = zeros(size(C));
C_offset = zeros(K, 1);
S = zeros(size(C));
A = full(A);
U = A'*Y;
V = A'*A;
aa = diag(V);   % squares of l2 norm all all components
sn =  zeros(1, K);
smin = zeros(1,K);
% kernel = obj.kernel;
kernel_pars = cell(K,1);
%% updating
ind_del = false(K, 1);
for miter=1:maxIter
    for k=1:K
        if  aa(k)==0
            C_raw(k, :) = C_raw(k, :)*0;
            C(k,:) = C(k, :)*0;
            S(k, :) = S(k, :)*0;
            ind_del(k) = true;
            continue;
        end
        if ind_del
            continue;
        end
        temp = C(k, :) + (U(k, :)-V(k, :)*C)/aa(k);
        %remove baseline and estimate noise level
        [b_hist, sn_hist] = estimate_baseline_noise(temp);
        b = mean(temp(temp<median(temp)));
        sn_psd = GetSn(temp);
        if sn_psd<sn_hist
            tmp_sn = sn_psd;
        else
            tmp_sn = sn_hist;
            b = b_hist;
        end

        temp = temp -b;
        sn(k) = tmp_sn;
        % deconvolution
        if obj.options.deconv_flag
            if strcmp(mode,'initiation')
                try
                    [ck, sk, deconv_options]= deconvolveCa(temp, deconv_options_0, 'sn', tmp_sn, 'maxIter', 2);
                    smin(k) = deconv_options.smin;
                    kernel_pars{k} = reshape(single(deconv_options.pars), 1, []);
                catch %if not deconvolved successfully
                    ck = max(0, temp);
                    sk=zeros(1,size(C,2));
                    if or(strcmp(deconv_options_0.type,'ar1'),strcmp(deconv_options_0.type,'kernel'))
                        kernel_pars{k}=single(0);
                    else
                        kernel_pars{k}=single(zeros(1,2));
                    end
                    ind_del(k) = true;
                end
            end
            
            if strcmp(mode,'massive')
                ck = temp;%max(0, temp);
                sk=zeros(1,size(C,2));
                if or(strcmp(deconv_options_0.type,'ar1'),strcmp(deconv_options_0.type,'kernel'))
                    kernel_pars{k}=single(0); 
                else
                    kernel_pars{k}=single(zeros(1,2));
                end
                ind_del(k) = true;
            end
            temp = temp - deconv_options.b; 
        else
            ck = max(0, temp);
        end                        
        C(k, :) = ck;   % save convolution kernels and deconvolution results (or not deconvolutioned)
        
        if sum(ck(2:end))==0
            ind_del(k) = true;
        end
        % save the spike count in the last iteration
        if miter==maxIter
            if obj.options.deconv_flag
                S(k, :) = sk;
            end
            C_raw(k, :) = temp;
            
        end
    end
end

obj.A = bsxfun(@times, A, sn);
obj.C = bsxfun(@times, C, 1./sn');
obj.C_raw = bsxfun(@times, C_raw, 1./sn');
obj.S = bsxfun(@times, S, 1./sn');
obj.P.kernel_pars =cell2mat(kernel_pars);
obj.P.smin = smin/sn;
obj.P.sn_neuron = sn;

if strcmp(mode,'initiation')
    if allow_deletion
        obj.delete(ind_del);
    end    
end

