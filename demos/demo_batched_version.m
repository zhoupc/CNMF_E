%% Batched CNMF-E version
% 1. (1) For each file, get A and C.
%    (2) Save in cell File. File(1).A=A;...
% 2. A=cat(2,File.A)
% 3. use partial extract_ac to get C
% 4. merge based on A C. We get A.
% 5. Use A to feed into each C same function as in 3.
% 6. Update A and C for each file. File(1).A=A;
% 7. Concatnate A and C, output neuron.A and neuron.C
%% 1
%%%% Input1 %%%%
datadir='/Volumes/data0-shared/elm/ProcessedCalciumData/sleep/7030/061117_5140F/';
kind='*CaELM*';
gSig=10;
gSiz=17;
min_pnr=6.5;
bg_neuron_ratio = 1;  % spatial range / diameter of neurons

Datadir=[datadir,kind];
filelist=dir(Datadir);

File(length(filelist)) = struct('HY',[],'A',[],'C',[],'C_raw',[],'S',[],'P',[],'options',[]);

parfor i= 1:length(filelist)
    nam=fullfile(datadir,filelist(i).name);
    demo_endoscope2
    File(i).HY=HY;
    File(i).A=neuron.A;
    File(i).options=neuron.options;
%     File(i).C=neuron.C;
%     File(i).C_raw=neuron.C_raw;
%     File(i).S=neuron.S;
%     File(i).P=neuron.P;
end

%% 2
A=cat(2,File.A);
Amask=A>0;
PNPN=sum(Amask,1); %pixelsNumberPerNeuron
Amask=reshape(Amask,[],1);
%% 3
parfor i= 1:length(filelist)
    HY=repmat(File(i).HY,size(A,2),1);    
    pInvolved = HY(Amask,:);
    neuron_p=mat2cell(pInvolved,PNPN,size(HY,2));    
    C_raw = cell2mat(cellfun(@mean, neuron_p,'UniformOutput',false))
    File(i).C_raw=C_raw; %%%%%%%%% Ybg
    
end
C_raw=cat(1,File.C_raw);
%% 4
[merged_ROIs, newIDs,A,~] = mergeAC(A,C_raw,merge_thr);
%% 5
Amask=A>0;
PNPN=sum(Amask,1); %pixelsNumberPerNeuron
Amask=reshape(Amask,[],1);
parfor i= 1:length(filelist)
    HY=repmat(File(i).HY,size(A,2),1);    
    pInvolved = HY(Amask,:);
    neuron_p=mat2cell(pInvolved,PNPN,size(HY,2));    
    C_raw = cell2mat(cellfun(@mean, neuron_p,'UniformOutput',false))
    File(i).C_raw=C_raw;    
    % denoise C_raw
    for j=1:size(C_raw,2)
        ci=C_raw(j,:);
        ai=A(:,j);
        [b, sn] = estimate_baseline_noise(ci);
        psd_sn = GetSn(ci);
        if sn>psd_sn
            sn =psd_sn;
            [ci, ~] = remove_baseline(ci, sn);
        else
            ci = ci - b;            
        end       
        ind_neg = (ci<-4*sn);
        ci(ind_neg) = rand(sum(ind_neg), 1)*sn;            
        % normalize the result
        ci = ci / sn;
        ai = ai * sn;
        % return results
        A(:,j)=ai;
        C_raw(j,:)=ci;
end








    
    
    
    
    
    
    