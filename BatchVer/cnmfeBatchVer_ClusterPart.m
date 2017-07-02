if running_on_cluster
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end

File(length(samplelist)) = struct('options',[],'Ysignal',[]); % pre-allocate for parfor loop. 
A0s=cell(1,length(samplelist));
%%% Running normal CNMF-E for each file
parfor i= 1:length(samplelist)
    Mode='initiation';
    picname=samplelist(i).name(namepattern) % For each file, save A's so you can roughly check what neuron is picked in which file. 
    name=fullfile(sampledir,samplelist(i).name);
    [A0s{i},File(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,FS,SSub,TSub,bg_neuron_ratio,name,Mode,picname,[],File(i),convolveType);
    fprintf('Sampling file number %.0f done\n', i);
end

%%% delete samples that have no neurons.
emptyA0s_ind=find(cellfun('isempty', A0s));
if ~isempty(emptyA0s_ind)
    warning(['sample file number ',num2str(emptyA0s_ind),' with name below has/have no neuron extracted in it.\n'])
    samplelist(emptyA0s_ind).name
    fprintf('Deleting these sample A0s.');
    samplelist_reduced=samplelist;
    samplelist_reduced(emptyA0s_ind)=[];
    A0s(emptyA0s_ind)=[];    
    %File(emptyA0s_ind)=[];
end
save([outputdir 'NormalsOFcnmfeBatchVer.mat'],'-v7.3')

%%% Order similar neurons in the same sequence in each file, not necessary,
%%% but nice to do. It is fast.
A0s=Over_Days_ResequenceA(A0s,correlation_thresh,max2max2nd,skewnessthresh);

A=cat(2,A0s{:}); Amask=A>0;  % This A is raw and plain, it is just all A's from all files concatnated together.
% Next, Use this A, in each file i, find C's corresponding to each A's found in file j.
ACS(length(samplelist)) = struct('Ain',[],'Cin',[],'STD',[]);
S_R=length(samplelist_reduced);
parfor i= 1:length(samplelist)
    Ain=[];
    Cin=[];
    STD=[];
    for j=1:S_R % parfor needs this
        Aj=A0s{j};
        ACS_temp=A2C2A(File(i), Aj, File(i).options);
        Ain = [Ain ACS_temp.Ain];
        Cin = [Cin; ACS_temp.Cin];
        STD=[STD ACS_temp.STD];
    end
    ACS(i).Ain=Ain;
    ACS(i).Cin=Cin;
    ACS(i).STD=STD;
end
save([outputdir 'PartOneOFcnmfeBatchVer.mat'],'-v7.3')
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