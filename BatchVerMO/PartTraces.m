function [PartC, AllC,kernel_pars, tmp_sn, boundary,neuron_batch,boundary_raw]=PartTraces(neuron_batch)
% This function uses one neurons's calcium traces from all time to estimate
% deconvolution parameters and apply it to all small traces extracted from seperate
% files. It concatenates files by noise-signal-noise so that the
% concatenated signals are coreherent.

% Modified by Shijie Gu, ShanghaiTech University; Fee Lab at BCS, MIT.
% Main dependence: deconvolveCa.mat in the original cnmf_e package.

K=size(neuron_batch(1).rawsignal,1);
PartC={}; PartC_raw={};
boundary={};
boundary_raw={};
kernel_pars={};

[AllC,~]=AllTraces(neuron_batch);
deconv_options_0 = neuron_batch(1).neuron.options.deconv_options;
tmp_sn=cell(1,K); b=cell(1,K);

for ni=1:K
    display(ni)
    C=[]; bound=[]; bound_raw=[];
    inpoint=max(AllC(ni,:),[],2)/50; % estimate noise level first point
    for fi=1:length(neuron_batch)
        C_single=neuron_batch(fi).rawsignal(ni,:);
        if range(C)/std(C)>3
            [b{ni}{fi}, tmp_sn{ni}{fi}] = estimate_baseline_noise(C_single);
        else
            b{ni}{fi} = mean(C_single(C_single<median(C_single(C_single~=0))));%%%%%%%%%%%
            tmp_sn{ni}{fi} = GetSn(C_single);
        end
        
        sig_temp_inpoint=find(C_single<inpoint,1,'first'); % estimate first noise level point
        sig_temp_outpoint=find(C_single<inpoint,1,'last'); % estimate last noise level point
        C_tmp=C_single(sig_temp_inpoint:sig_temp_outpoint);
        C_tmp=(C_tmp-b{ni}{fi})./tmp_sn{ni}{fi};
        C=[C C_tmp];                 % concatenate temporal signals
        bound=[bound size(C_tmp,2)]; %for debug
    end

    %% modification from updateTemporal
%     if range(C)/std(C)>6
%         [b, tmp_sn] = estimate_baseline_noise(C);
%     else
%        b = mean(C(C<median(C(C~=0))));%%%%%%%%%%%
%        tmp_sn = GetSn(C);
%     end
%         end
%     b = mean(C(C<median(C(C~=0))));%%%%%%%%%%%
%     tmp_sn = GetSn(C);
    %C(C<0)=0;
    PartC_raw{ni}=C;
    tmp_sn_all=mean(cat(2,tmp_sn{ni}{:}));
    [~, ~, deconv_options]= deconvolveCa(C, deconv_options_0, 'sn', 1, 'maxIter', 2);
    kernel_pars{ni} = reshape(single(deconv_options.pars), 1, []);
    boundary{ni}=cumsum(bound); %for debug
    
    
    tmp=[];
    for fi=1:length(neuron_batch)        
        [sig_tmp, ~, ~]= deconvolveCa(neuron_batch(fi).rawsignal(ni,:), deconv_options_0,'pars', (kernel_pars{ni})','sn', tmp_sn{ni}{fi}, 'maxIter', 2);
        neuron_batch(fi).signal(ni,:)=sig_tmp./tmp_sn{ni}{fi};
        tmp=[tmp neuron_batch(fi).signal(ni,:)];
        bound_raw=[bound_raw size(neuron_batch(fi).signal(ni,:),2)];
    end
    PartC{ni}=tmp;
    boundary_raw{ni}=cumsum(bound_raw);
end    
    
end