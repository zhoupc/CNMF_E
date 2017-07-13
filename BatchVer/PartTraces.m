function [PartC, kernel_pars, sn, boundary,neuron_batch]=PartTraces(neuron_batch)
K=size(neuron_batch(1).signal,1);
PartC={};
boundary={};
kernel_pars={};

[AllC,~]=AllTraces(neuron_batch);
inpoint=max(AllC(1,:),[],2)/50;
deconv_options_0 = neuron_batch(1).neuron.options.deconv_options;

for ni=1:K
    C=[]; bound=[];
    for fi=1:length(neuron_batch)
        sig_temp_inpoint=find(neuron_batch(fi).signal(ni,:)<inpoint,1,'first');
        sig_temp_outpoint=find(neuron_batch(fi).signal(ni,:)<inpoint,1,'last');
        C_tmp=neuron_batch(fi).signal(ni,sig_temp_inpoint:sig_temp_outpoint);
        C=[C C_tmp];
        bound=[bound size(C_tmp,2)]; %for debug
    end

    %% modification from updateTemporal
%             if range(temp)/std(temp)>6
%             [b, tmp_sn] = estimate_baseline_noise(temp);
%         else
%             b = mean(temp(temp<median(temp)));
%             tmp_sn = GetSn(temp);
%         end
    b = mean(C(C<median(C(C~=0))));%%%%%%%%%%%
    tmp_sn = GetSn(C);
    C = C -b;
    %C(C<0)=0;
    sn{ni} = tmp_sn;
    [PartC{ni}, ~, deconv_options]= deconvolveCa(C, deconv_options_0, 'sn', tmp_sn, 'maxIter', 2);
    kernel_pars{ni} = reshape(single(deconv_options.pars), 1, []);
    boundary{ni}=cumsum(bound); %for debug
    
    for fi=1:length(neuron_batch)
        [neuron_batch(fi).signal(ni,:), ~, ~]= deconvolveCa(neuron_batch(fi).signal(ni,:), deconv_options_0,'pars', (kernel_pars{ni})','sn', sn{ni}, 'maxIter', 2);
    end
end    
    
end