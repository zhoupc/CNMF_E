function C_ = deconvTemporal(obj, use_parallel)
if ~exist('use_parallel', 'var')||isempty(use_parallel)
    use_parallel = true;
end
C_raw_ = obj.C_raw;
K = size(C_raw_, 1);
C_raw_ = mat2cell(C_raw_, ones(K,1), size(C_raw_,2));
C_ = cell(K,1);
S_ = C_;
kernel_pars = cell(K, 1);
sn = cell(K,1);
num_per_row = 100;
for m=1:K
    fprintf('|');
    if mod(m, num_per_row)==0
        fprintf('\n');
    end
end
fprintf('\n');
deconv_options = obj.options.deconv_options;
if use_parallel
    tmp_flag = false(K,1); 
    ind = randi(K, ceil(K/num_per_row), 1); 
    tmp_flag(ind) = true; 
    tmp_flag = num2cell(tmp_flag); 
    parfor k=1:size(C_raw_,1)
        ck_raw = C_raw_{k};
        
        %remove baseline and estimate noise level. estiamte the noise
        %using two methods: psd and histogram. choose the one with
        %smaller value.
        [b_hist, sn_hist] = estimate_baseline_noise(ck_raw);

        baseline_ = mean(ck_raw(ck_raw<median(ck_raw)));
        sn_psd = GetSn(ck_raw);
        if sn_psd<sn_hist
            tmp_sn = sn_psd;
        else
            tmp_sn = sn_hist;
            baseline_ = b_hist;
        end
        
        % subtract the baseline
        ck_raw = ck_raw -baseline_;
        sn{k} = tmp_sn;
        
        % deconvolution
        [ck, sk, tmp_options]= deconvolveCa(ck_raw, deconv_options, 'maxIter', 2, 'sn', tmp_sn);
        if sum(abs(ck))==0
            ck = ck_raw;
        end
        C_{k} = reshape(ck, 1, []);
        S_{k} = reshape(sk, 1, []);
        kernel_pars{k} = reshape(tmp_options.pars, 1, []);
        C_raw_{k} = ck_raw - tmp_options.b;
        
        fprintf('.');
        if tmp_flag{k}
            fprintf('\n');
        end
    end
else
    for k=1:size(C_raw_,1)
        ck_raw = C_raw_{k};
        
        %remove baseline and estimate noise level. estiamte the noise
        %using two methods: psd and histogram. choose the one with
        %smaller value.
        [b_hist, sn_hist] = estimate_baseline_noise(ck_raw);
        baseline_ = mean(ck_raw(ck_raw<median(ck_raw)));
        sn_psd = GetSn(ck_raw);
        if sn_psd<sn_hist
            tmp_sn = sn_psd;
        else
            tmp_sn = sn_hist;
            baseline_ = b_hist;
        end
        
        % subtract the baseline
        ck_raw = ck_raw -baseline_;
        sn{k} = tmp_sn;
        
        % deconvolution
        [ck, sk, tmp_options]= deconvolveCa(ck_raw, deconv_options, 'maxIter', 2, 'sn', tmp_sn);
        
        if sum(abs(ck))==0
            ck = ck_raw;
        end
        C_{k} = reshape(ck, 1, []);
        S_{k} = reshape(sk, 1, []);
        kernel_pars{k} = reshape(tmp_options.pars, 1, []);
        C_raw_{k} = ck_raw - tmp_options.b;
        fprintf('.'); 
        if mod(k,num_per_row)==0
            fprintf('\n');
        end
    end
end
fprintf('\n');
C_ = cell2mat(C_);
obj.C = C_;
obj.S = cell2mat(S_);
obj.P.kernel_pars = cell2mat(kernel_pars);
obj.P.neuron_sn = cell2mat(sn);
end
