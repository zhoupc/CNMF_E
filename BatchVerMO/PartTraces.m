function [PartC, AllC,kernel_pars, tmp_sn, boundary,neuron_batch,boundary_raw]=PartTraces(neuron_batch)
% This function uses one neurons's calcium traces the temporal traces from all files to estimate
% deconvolution parameters and apply it to all small traces extracted from seperate
% files. It concatenates files by noise-signal-noise so that the
% concatenated signals are coreherent.

% Outputs: 
%   PartC: all the neurons' temporal traces concatenated, 
%       neuron_batch.signal. neurons are put in rows.
%   AllC: all the neurons' temporal traces concatenated. 
%       Similar to PartC, but no deconvolved, only rescaled(normalized by noise)
%       Neurons are put in rows.
%   kernel_pars: cell array, length of neuron number. Each cell has the
%       parameters, 1*N parameter, based on your choice of deconvolution.
%   tmp_sn: cell array, length of neuron number. Each cell is a cell array
%       of length of file number. Cell{i}{j} is the noise level of neuron
%       {i} in file {j}
%   PartC_raw(omited, but can always get): 
%       Cell array, of length of neuron number. Each cell contains the noise-signal-noise
%       concatenated temporal traces.
%   boundary: cell array, of length of neuron number. Each cell is the time
%       boundary for PartC_raw.
%   neuron_batch: field 'signal' is filled in.
%   boundary_raw: similar to boundary, although it is the boundary for
%       PartC.

%   
% Modified by Shijie Gu, ShanghaiTech University; Fee Lab at BCS, MIT.
% Main dependence: deconvolveCa.mat in the original cnmf_e package.

K=size(neuron_batch(1).rawsignal,1);
PartC=[]; PartC_raw={};
boundary={};
boundary_raw=[];
kernel_pars={};

deconv_options_0 = neuron_batch(1).neuron.options.deconv_options;
tmp_sn=cell(1,K); b=cell(1,K);

for ni=1:K
    display(['neuron' ni])
    C=[]; bound=[]; bound_raw=[];
    
    for fi=1:length(neuron_batch)
        C_single=neuron_batch(fi).rawsignal(ni,:);
        if range(C)/std(C)>3
            [b{ni}{fi}, tmp_sn{ni}{fi}] = estimate_baseline_noise(C_single);
        else
            b{ni}{fi} = mean(C_single(C_single<median(C_single(C_single~=0))));%%%%%%%%%%%
            tmp_sn{ni}{fi} = GetSn(C_single);
        end
        C_single=(C_single-b{ni}{fi})./tmp_sn{ni}{fi};  
        neuron_batch(fi).signal(ni,:)=C_single;
    end
    [AllC,~]=AllTraces(neuron_batch,'signal');
    inpoint=max(AllC(ni,:),[],2)/50; % estimate noise level first point
    for fi=1:length(neuron_batch) 
        C_single=neuron_batch(fi).signal(ni,:);
        sig_temp_inpoint=find(C_single<inpoint,1,'first'); % estimate first noise level point
        sig_temp_outpoint=find(C_single<inpoint,1,'last'); % estimate last noise level point
        C_tmp=C_single(sig_temp_inpoint:sig_temp_outpoint);
        C=[C C_tmp];                 % concatenate temporal signals
        bound=[bound size(C_tmp,2)]; %for debug
    end

    PartC_raw{ni}=C;            %for debug
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
    PartC(ni,:)=tmp;
    boundary_raw=cumsum(bound_raw);
end    
    
end