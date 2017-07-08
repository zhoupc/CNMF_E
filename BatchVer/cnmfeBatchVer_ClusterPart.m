%% C Run All below on cluster!
%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end

%% 1. Run normal CNMF-E for each file
File(length(samplelist)) = struct('options',[],'Ysignal',[]); % pre-allocate for parfor loop. 
A0s=cell(1,length(samplelist));
parfor i= 1:length(samplelist)
    Mode='initiation';
    picname=samplelist(i).name(namepattern) % For each file, save A's so you can roughly check what neuron is picked in which file. 
    name=fullfile(sampledir,samplelist(i).name);
    [A0s{i},File(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,min_pixel,bd,FS,SSub,TSub,bg_neuron_ratio,name,Mode,picname,[],File(i),convolveType,merge_thr);
    fprintf('Sampling file number %.0f done\n', i);
end

%%% delete samples that have no neurons.
emptyA0s_ind=find(cellfun('isempty', A0s));
if ~isempty(emptyA0s_ind)
    warning(['sample file number ',num2str(emptyA0s_ind),' with name below has/have no neuron extracted in it.\n'])
    samplelist(emptyA0s_ind).name
    fprintf('Deleting these sample A0s.');

    samplelist(emptyA0s_ind)=[];
    A0s(emptyA0s_ind)=[];    
    File(emptyA0s_ind)=[];
end
save([outputdir 'NormalsOFcnmfeBatchVer.mat'],'-v7.3')

%%% Order similar neurons in the same sequence in each file, not necessary,
%%% but nice to do. It is fast.
A0s=Over_Days_ResequenceA(A0s,correlation_thresh,max2max2nd,skewnessthresh);

%% 2. Next, Use this A, in each file i, find C's corresponding to each A's found in file j.
ACS(length(samplelist)) = struct('Ain',[],'Cin',[],'STD',[]);
S_R=length(samplelist);
parfor i= 1:S_R
    Ain=[]; Cin=[]; STD=[];
    for j=1:S_R % parfor needs this
        Aj=A0s{j};
        ACS_temp=A2C2A(File(i), Aj, File(i).options);
        Ain = [Ain ACS_temp.Ain]; Cin = [Cin; ACS_temp.Cin]; STD=[STD ACS_temp.STD];
    end
    ACS(i).Ain=Ain; ACS(i).Cin=Cin; ACS(i).STD=STD;
end
save([outputdir 'PartOneOFcnmfeBatchVer.mat'],'-v7.3')
%% 3 Merge similar neurons
%%% Merge similar neurons based on spatial AND temporal correlation
% C_all=cat(2,ACS.Cin);
% Amask_temp=cat(2,A0s{1:2});
% Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3)); %only use central part for merging.

% C_temp=C_all(1:size(Amask_temp,2),:);
% [Amask_temp,C_temp,ACS] = mergeAC(Amask_temp,C_temp,ACS,merge_thr_2);
% 
% merge2start=1+size(A0s{1},2)+size(A0s{2},2);
% for i=3:length(samplelist)
%     Amask_temp=cat(2,Amask_temp,A0s{i});
%     Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3));
%     
%     C_temp=[C_temp;C_all(merge2start:merge2start+size(A0s{i},2)-1,:)];
% 
%     [Amask_temp,C_temp,ACS] = mergeAC(Amask_temp,C_temp,ACS,merge_thr_2);
%     merge2start=merge2start+size(A0s{i},2);
% end

Amask_temp=cat(2,A0s{:});
Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3)); %only use central part for merging.
%C_all=cat(2,ACS.Cin);

[Afinal,MC,newIDs,merged_ROIs] = mergeAC(Amask_temp,ACS,merge_thr_2);
% [size1,~]=cellfun(@size,newIDs);
% ACS(size1~=1)
save([outputdir 'commonAcnmfeBatchVer.mat'],'-v7.3')
%% 4 Collapsing A's
% As=cell(1,length(samplelist));
% STDs=cell(1,length(samplelist));
% parfor i=1:length(samplelist)
%     As{i}=ACS(i).Ain;
%     STDs{i}=ACS(i).STD;    
% end
% Afinal=ReducingA(As,STDs); % the format for cell input is designed for potential other applications.
% save([outputdir 'uniqueAcnmfeBatchVer.mat'],'-v7.3')
%% 4.5 Determine Afinal that will be used to extract C's in each file.

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
newIDs=newIDs(nz_ind);

save([outputdir 'AfinalcnmfeBatchVer.mat'],'-v7.3')

%% 5 "massive" procedure: Extract A from each file
neuron_batch(length(filelist)) = struct('ind_del',[],'signal',[],'FileOrigin',[],'neuron',[]);

parfor i= 1:length(filelist)  
    mode='massive';
    nam=fullfile(datadir,filelist(i).name);    
    [~,neuron_batch(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,min_pixel,bd,FS,SSub,TSub,bg_neuron_ratio,nam,mode,[],Afinal,neuron_batch(i),convolveType,merge_thr);
    neuron_batch(i).FileOrigin=filelist(i); % save origin(filelist)
end
fprintf('Massive extraction done.');
save([outputdir 'MassivecnmfeBatchVer.mat'],'-v7.3')

%% 6 Partition between those neurons found in each file and those not. Save results.
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