function [PartC, kernel_pars, sn, boundary,neuron_batch]=PartTraces(neuron_batch)
% This function uses one neurons's calcium traces from all time to estimate
% deconvolution parameters and apply it to all small traces extracted from seperate
% files. It concatenates files by noise-signal-noise so that the
% concatenated signals are coreherent.

% Modified by Shijie Gu, ShanghaiTech University; Fee Lab at BCS, MIT.
% Main dependence: deconvolveCa.mat in the original cnmf_e package.

K=size(neuron_batch(1).rawsignal,1);
PartC={};
boundary={};
kernel_pars={};

[AllC,~]=AllTraces(neuron_batch);
inpoint=max(AllC(1,:),[],2)/50; % estimate noise level first point
deconv_options_0 = neuron_batch(1).neuron.options.deconv_options;

for ni=1:K
    display(ni)
    C=[]; bound=[];
    for fi=1:length(neuron_batch)
        sig_temp_inpoint=find(neuron_batch(fi).rawsignal(ni,:)<inpoint,1,'first'); % estimate first noise level point
        sig_temp_outpoint=find(neuron_batch(fi).rawsignal(ni,:)<inpoint,1,'last'); % estimate last noise level point
        C_tmp=neuron_batch(fi).rawsignal(ni,sig_temp_inpoint:sig_temp_outpoint);
        C=[C C_tmp];                 % concatenate temporal signals
        bound=[bound size(C_tmp,2)]; %for debug
    end

    %% modification from updateTemporal
    if range(C)/std(C)>6
        [b, tmp_sn] = estimate_baseline_noise(C);
    else
        b = mean(C(C<median(C(C~=0))));%%%%%%%%%%%
        tmp_sn = GetSn(C);
    end
%         end
%     b = mean(C(C<median(C(C~=0))));%%%%%%%%%%%
%     tmp_sn = GetSn(C);
    C = C -b;
    %C(C<0)=0;
    sn{ni} = tmp_sn;
    [PartC{ni}, ~, deconv_options]= deconvolveCa(C, deconv_options_0, 'sn', tmp_sn, 'maxIter', 2);
    kernel_pars{ni} = reshape(single(deconv_options.pars), 1, []);
    boundary{ni}=cumsum(bound); %for debug
    
    for fi=1:length(neuron_batch)
        [neuron_batch(fi).signal(ni,:), ~, ~]= deconvolveCa(neuron_batch(fi).rawsignal(ni,:), deconv_options_0,'pars', (kernel_pars{ni})','sn', sn{ni}, 'maxIter', 2);
    end
end    
    
end