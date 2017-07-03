%% 2 Merge similar neurons
%%% Merge similar neurons based on spatial AND temporal correlation
C_all=cat(2,ACS.Cin);

Amask_temp=cat(2,A0s{1:2})>0;
C_temp=C_all(1:size(Amask_temp,2),:);
[Amask_temp,C_temp,ACS] = mergeAC(Amask_temp,C_temp,ACS,merge_thr);

merge2start=1+size(A0s{1},2)+size(A0s{2},2);
for i=3:length(samplelist_reduced)
    Amask_temp=cat(2,Amask_temp,A0s{i})>0;
    
    C_temp=[C_temp;C_all(merge2start:merge2start+size(A0s{i},2)-1,:)];

    [Amask_temp,C_temp,ACS] = mergeAC(Amask_temp,C_temp,ACS,merge_thr);
    merge2start=merge2start+size(A0s{i},2);
end

save([outputdir 'commonAcnmfeBatchVer.mat'],'-v7.3')
%% 3 Determine uniqueA's
As=cell(1,length(samplelist));
STDs=cell(1,length(samplelist));
parfor i=1:length(samplelist)
    As{i}=ACS(i).Ain;
    STDs{i}=ACS(i).STD;    
end
Afinal=ReducingA(As,STDs);
save([outputdir 'uniqueAcnmfeBatchVer.mat'],'-v7.3')
%% 5 Determine Afinal that will be used to extract C's in each file.

%%% Some processes making Afinal nicer, modified from Pengcheng Zhou's
%%% idea.
for i=1:size(Afinal,2)
    ai=Afinal(:,i);
    temp = full(ai>quantile(ai, 0.5, 1));
    ai(~temp(:)) = 0;
    Afinal(:,i)=ai;
end

% Just in case some all zero A's got passed to this stage.
nz_ind=any(Afinal);
Afinal=Afinal(:,nz_ind);
save([outputdir 'AfinalcnmfeBatchVer.mat'],'-v7.3')
%% 6 "massive" procedure: Extract A from each file
neuron_batch(length(filelist)) = struct('ind_del',[],'signal',[],'FileOrigin',[],'neuron',[]);

parfor i= 1:length(filelist)  
    mode='massive';
    nam=fullfile(datadir,filelist(i).name);    
    [~,neuron_batch(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,FS,SSub,TSub,bg_neuron_ratio,nam,mode,[],Afinal,neuron_batch(i),convolveType);
    neuron_batch(i).FileOrigin=filelist(i); % save origin(filelist)
end
fprintf('Massive extraction done.');
save([outputdir 'MassivecnmfeBatchVer.mat'],'-v7.3')

%neuron(length(filelist)) = struct('signal',[],'filelist',[]);
%%% Partition between those neurons found in each file and those not.
ind_del_final_cat=cat(2,neuron_batch.ind_del);
ind_del_final=any(ind_del_final_cat,2);
parfor i= 1:length(filelist)
    neuron_batch(i).neuron.A=[neuron_batch(i).neuron.A(:,~ind_del_final) neuron_batch(i).neuron.A(:,ind_del_final)];
    fprintf('A extraction done\n');
    neuron_batch(i).neuron.C=[neuron_batch(i).neuron.C(~ind_del_final,:);neuron_batch(i).neuron.C(ind_del_final,:)];
    neuron_batch(i).neuron.C_raw=[neuron_batch(i).neuron.C_raw(~ind_del_final,:);neuron_batch(i).neuron.C_raw(ind_del_final,:)];
    fprintf('C extraction done\n');
            %%% save each data i's signal into neuron(i).signal and
    for j=1:size(neuron_batch(i).neuron.A,2)
        jA=neuron_batch(i).neuron.A(:,j);
        jC=neuron_batch(i).neuron.C(j,:);
        neuron_batch(i).signal(j,:)=median(jA(jA>0)*jC);
    end        
    fprintf('neuron_batch %.0f extraction done\n', i);
end
fprintf('First %.0f neurons are successfully deconvolved in each files while those after that are missing in some files\n', sum(~ind_del_final));
fprintf('ALL extractions done.\n');
eval(sprintf('save %sCNMFE_BatchVer.mat %s -v7.3', outputdir, 'neuron_batch'));
fprintf('ALL data saved, check them out!');