%% C Run All below on cluster!
if strcmp(Version,'MoBatchVer')
    fprintf('Running BatchVer with Motion Correction Later.')
elseif strcmp(Version,'BatchVer')
    fprintf('Simple BatchVer')    
end
if ~exist(outputdirDetails,'dir')
    mkdir(outputdirDetails)
end
cd(outputdirDetails)

%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end

%% 1. Run normal CNMF-E for each file
File(length(samplelist)) = struct('options',[],'Ysignal',[],'neuron',[],'Ybg',[]); % Ysignal and Ybg for video making.
                                                                                   % pre-allocate for parfor loop. 
A0s=cell(1,length(samplelist));
parfor i= 1:length(samplelist)
    Mode='initiation';
    picname=samplelist(i).name %(namepattern) % For each file, save A's so you can roughly check what neuron is picked in which file. 
    name=fullfile(sampledir,samplelist(i).name);
    [A0s{i},File(i)]=demo_endoscope2(bg_neuron_ratio,merge_thr,with_dendrites,K,sframe,num2read,...
                                   name,neuron_full,Mode,picname,File(i),[],...
                                   thresh_detecting_frames);                               
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

%%% Order similar neurons in the same sequence in each file, not necessary,
%%% but nice to do. It is fast.
A0s=Over_Days_ResequenceA(A0s,correlation_thresh,max2max2nd,skewnessthresh);
save([outputdirDetails 'EachFilecnmfeBatchVer.mat'],'-v7.3')

if strcmp(Version,'MoBatchVer')
    eval(sprintf('save %s%0.f_cnmfe_BatchVer_PartI_File.mat %s -v7.3', outputdir, daynum, 'File'));
    fprintf('Variable ''File'', which includes fields of ''options'' ''Ysignal'' ''neuron'' ''Ybg'' saved.');
elseif strcmp(Version,'BatchVer')
    Vars = {'File'}; Vars=Vars';
    eval(sprintf('save %s%0.f_cnmfe_BatchVer_ClusterPartI.mat %s -v7.3', outputdir, daynum, strjoin(Vars)));
end
File = rmfield(File,{'Ybg','neuron'});

%% 2. Next, Use this A, in each file i, find C's corresponding to each A's found in file j.
ACS(length(samplelist)) = struct('Cin',[],'Cin_raw',[],'STD',[]);
S_R=length(samplelist);
parfor i= 1:S_R
    Cin=[]; Cin_raw=[]; STD=[];
    for j=1:S_R % parfor needs S_R rather than length(samplelist)
        Aj=A0s{j};
        ACS_temp=A2C(File(i).Ysignal, Aj, File(i).options);
        Cin = [Cin; ACS_temp.Cin]; STD=[STD ACS_temp.STD]; Cin_raw=[Cin_raw; ACS_temp.Cin_raw];
    end
    ACS(i).Cin=Cin; ACS(i).Cin_raw=Cin_raw; ACS(i).STD=STD;
end
save([outputdirDetails 'ACScnmfeBatchVer.mat'],'-v7.3')

%% 3 Merge similar neurons

A_temp=cat(2,A0s{:});
%Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3)); %only use central part for merging.
[Afinal,MC,newIDs,merged_ROIs,close_ind,real_ind] = mergeAC(A_temp,ACS,merge_thr_2,dmin,File(1).options.d1,File(1).options.d2);

if strcmp(Version,'MoBatchVer'); save([outputdirDetails 'commonAcnmfeBatchVer.mat'],'-v7.3')
elseif strcmp(Version,'BatchVer');save([outputdir 'commonAcnmfeBatchVer.mat'],'-v7.3')
end

%% 4 Determine Afinal that will be used to extract C's in each file.

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


Apicname=sprintf('%.0fAfinal',daynum);
if strcmp(Version,'MoBatchVer')
    ColorAllNeurons(Afinal,File(1).options.d1,File(1).options.d2,[Apicname ' PNR=' num2str(neuron_full.options.min_pnr)],outputdirDetails);
    Vars = {'Afinal';'samplelist'}; Vars=Vars';
    eval(sprintf('save %s%0.f_cnmfe_BatchVer_PartI_Afinalsam.mat %s -v7.3', outputdir, daynum, strjoin(Vars)));
    fprintf('cnmfe_BatchVer_for motion Part1 data saved, check them out!');
    return
elseif strcmp(Version,'BatchVer')
    ColorAllNeurons(Afinal,File(1).options.d1,File(1).options.d2,Apicname,outputdir);
    Vars = {'Afinal';'samplelist'}; Vars=Vars';
    eval(sprintf('save %s%0.f_cnmfe_BatchVer_ClusterPartI.mat %s -append', outputdir, daynum, strjoin(Vars)));
end

% The following will be executed for cnmf_e(BatchVer), without motion
% correction.
%% 5 "massive" procedure: Extract A from each file
neuron_batch(length(samplelist)) = struct('ind_del',[],'rawsignal',[],'signal',[],'DeconvSpiketrain',[],'FileOrigin',[],'neuron',[],'C',[],'C_raw',[]);

parfor i= 1:length(samplelist)  
    mode='massive';
    nam=cell(1,2);
    nam{1}=fullfile(datadir,samplelist(i).name);
    nam{2}=File(i).Ysignal;
    [~,neuron_batch(i)]=demo_endoscope2(bg_neuron_ratio,[],with_dendrites,K,sframe,num2read,...
                                   nam,neuron_full,mode,[],neuron_batch(i),Afinal,...
                                   thresh_detecting_frames);
    neuron_batch(i).FileOrigin=filelist(i); % save origin(filelist)
end
neuron_batch = rmfield(neuron_batch,{'C','C_raw'});
fprintf('Massive extraction in each file done.');

%% 6 Save A*C

parfor i= 1:length(samplelist)
    for j=1:size(neuron_batch(i).neuron.A,2)
        jA=neuron_batch(i).neuron.A(:,j);
        jC=neuron_batch(i).neuron.C_raw(j,:);
        neuron_batch(i).rawsignal(j,:)=median(jA(jA>0)*jC);
    end        
    fprintf('neuron_batch %.0f extraction done\n', i);
end

%% 7 deconvolve signal
[~, ~, ~, ~,~,neuron_batch,boundary_raw]=PartTraces(neuron_batch);
Vars = {'neuron_batch';'boundary_raw'}; Vars=Vars';
fprintf('ALL extractions done.\n');
eval(sprintf('save %sCNMFE_BatchVer.mat %s -v7.3', outputdir, strjoin(Vars)));
fprintf('ALL data saved, check them out!');
