if running_on_cluster
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end

File(length(samplelist)) = struct('options',[],'Y',[],'Ysignal',[]); % pre-allocate for parfor loop. 
A0s=cell(1,length(samplelist));
%%% Running normal CNMF-E for each file
parfor i= 1:length(samplelist)
    Mode='initiation';
    picname=samplelist(i).name(namepattern) % For each file, save A's so you can roughly check what neuron is picked in which file. 
    name=fullfile(sampledir,samplelist(i).name);
    [A0s{i},File(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,FS,SSub,TSub,bg_neuron_ratio,name,Mode,picname,[],File(i),convolveType);
    fprintf('Sampling file number %.0f done\n', i);
end

save([outputdir 'NormalsOFcnmfeBatchVer.mat'])

%%% Order similar neurons in the same sequence in each file, not necessary,
%%% but nice to do. It is fast.
[ns_storage_1]=Over_Days_findAnn(A0s,correlation_thresh,max2max2nd,skewnessthresh);
for i= 1:length(samplelist)
    A0=A0s{i};
    corr_ind=ns_storage_1(:,i);
    unique_ind=setdiff(1:size(A0,2),corr_ind);
    A0s{i}=[A0(:,corr_ind) A0(:,unique_ind)];
end

A=cat(2,A0s{:}); Amask=A>0;  % This A is raw and plain, it is just all A's from all files concatnated together.
% Next, Use this A, in each file i, find C's corresponding to each A's found in file j.
ACS(length(samplelist)) = struct('Ain',[],'Cin',[],'STD',[]);
parfor i= 1:length(samplelist)
    for j=1:length(samplelist)
        Aj=A0s{j};
        ACS(i)=A2C2A(ACS(i), File(i), Aj, File(i).options);
    end
end
save([outputdir 'PartOneOFcnmfeBatchVer.mat'])
%% 2 Determine commonA's
%%% Merge similar neurons based on spatial AND temporal correlation
[~,commonA,commonC,ind_del] = mergeAC(Amask,ACS,merge_thr);
save([outputdir 'commonAcnmfeBatchVer.mat'],'-v7.3')
%% 3 Determine uniqueA's
Aunique=cell(1,length(samplelist));
STDunique=cell(1,length(samplelist));
parfor i=1:length(samplelist)
    FileA=ACS(i).Ain;
    FileSTD=ACS(i).STD;
    Aunique{i}=FileA(:,~ind_del);
    STDunique{i}=FileSTD(~ind_del);    
end
weightedA=ReducingA(Aunique,STDunique);
save([outputdir 'uniqueAcnmfeBatchVer.mat'],'-v7.3')
%% 5 Determine Afinal that will be used to extract C's in each file.
Afinal=cat(2,commonA,weightedA);

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
FILE(length(filelist)) = struct('A',[],'C',[],'ind_del',[]);

parfor i= 1:length(filelist)  
    mode='massive';
    nam=fullfile(datadir,filelist(i).name);    
    [~,FILE(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,FS,SSub,TSub,bg_neuron_ratio,nam,mode,[],Afinal,FILE(i),convolveType);    
end
fprintf('Massive extraction done.');
save([outputdir 'MassivecnmfeBatchVer.mat'],'-v7.3')

neuron(length(filelist)) = struct('signal',[],'filelist',[]);
%%% Partition between those neurons found in each file and those not.
ind_del_final_cat=cat(2,FILE.ind_del);
ind_del_final=any(ind_del_final_cat,2);
parfor i= 1:length(filelist)
    nA=FILE(i).A;
    FILE(i).A=[nA(:,~ind_del_final) nA(:,ind_del_final)];
    fprintf('A extraction done');
    nC=FILE(i).C;
    FILE(i).C=[nC(~ind_del_final,:);nC(ind_del_final,:)];
    fprintf('C extraction done');
    %%% save each data i's signal and origin(filelist) into neuron(i).signal and
    %%% neuron(i).filelist
    for j=1:size(FILE(i).A,2)
        jA=FILE(i).A(:,j);
        jC=FILE(i).C(j,:);
        neuron(i).signal(j,:)=median(jA(jA>0)*jC);
        neuron(i).filelist=filelist(i);
    end        
    fprintf('FILE and "neuron" %.0f extraction done\n', i);
end
fprintf('First %.0f neurons are found in each files while those after that are missing in some files', sum(~ind_del_final));
fprintf('ALL extractions done');
eval(sprintf('save %sCNMFE_BatchVer.mat %s -v7.3', outputdir, 'neuron FILE'));
fprintf('ALL data saved, check them out!');