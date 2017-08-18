%% C Run All below on cluster!
if strcmp(Version,'MoBatchVer')
    fprintf('Running BatchVer For Motion')
    if ~exist(outputdirDetails,'dir')
        mkdir(outputdirDetails)
    end
    cd(outputdirDetails)
elseif strcmp(Version,'BatchVer')
    fprintf('Simple BatchVer')
end

%% 0. Get cluster ready
if running_on_cluster % some procedures making cluster use robust
    [~, ~, ~] = maybe_spawn_workers(workersnum); 
    init_par_rng(2016);
end

%% 1. Run normal CNMF-E for each file

load([outputdirDetails 'ACScnmfeBatchVer.mat'])
merge_thr_2=[0.6,0.65];
% %% 2. Next, Use this A, in each file i, find C's corresponding to each A's found in file j.
% ACS(length(samplelist)) = struct('Cin',[],'Cin_raw',[],'STD',[]);
% S_R=length(samplelist);
% parfor i= 1:S_R
%     Cin=[]; Cin_raw=[]; STD=[];
%     for j=1:S_R % parfor needs this
%         Aj=A0s{j};
%         ACS_temp=A2C2A(File(i), Aj, File(i).options);
%         Cin = [Cin; ACS_temp.Cin]; STD=[STD ACS_temp.STD]; Cin_raw=[Cin_raw; ACS_temp.Cin_raw];
%     end
%     ACS(i).Cin=Cin; ACS(i).Cin_raw=Cin_raw; ACS(i).STD=STD;
% end
% % outputdir_video='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/videos/';
% % MakingVideos(File,File(1).options.d1,File(1).options.d2,num2str(daynum),outputdir_video)
% save([outputdirDetails 'ACScnmfeBatchVer.mat'],'-v7.3')
%% 3 Merge similar neurons

Amask_temp=cat(2,A0s{:});
%Amask_temp=bsxfun(@gt,Amask_temp,quantile(Amask_temp,0.3)); %only use central part for merging.
[Afinal,MC,newIDs,merged_ROIs,close_ind,real_ind] = mergeAC(Amask_temp,ACS,merge_thr_2,5,File(1).options.d1,File(1).options.d2);

if strcmp(Version,'MoBatchVer')
    save([outputdirDetails 'commonAcnmfeBatchVer.mat'],'-v7.3')
elseif strcmp(Version,'BatchVer')
    save([outputdir 'commonAcnmfeBatchVer.mat'],'-v7.3')
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
    ColorAllNeurons(Afinal,File(1).options.d1,File(1).options.d2,Apicname,outputdirDetails);
    Vars = {'Afinal';'samplelist'}; Vars=Vars';
    eval(sprintf('save %s%0.f_cnmfe_BatchVer_PartI_Afinalsam.mat %s -v7.3', outputdir, daynum, strjoin(Vars)));
    eval(sprintf('save %s%0.f_cnmfe_BatchVer_PartI_File.mat %s -v7.3', outputdir, daynum, 'File'));
    fprintf('cnmfe_BatchVer_for motion Part1 data saved, check them out!');

        outputdir_video='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/videos/';
        d1=File(1).options.d1; d2=File(1).options.d2;
        %MakingVideos(File,d1,d2,num2str(daynum),outputdir_video)
        clear ACS
        %MakingVideos([],d1,d2,num2str(daynum),outputdir_video,true,datadir,filelist)
        MakingVideos(File,d1,d2,num2str(daynum),outputdir,datadir,filelist,1)
        fprintf('Videos saved, check them out!');

    return
elseif strcmp(Version,'BatchVer')
    ColorAllNeurons(Afinal,File(1).options.d1,File(1).options.d2,Apicname,outputdir);
    Vars = {'Afinal';'samplelist';'File'}; Vars=Vars';
    eval(sprintf('save %s%0.f_cnmfe_BatchVer_ClusterPartI.mat %s -v7.3', Aoutputdir, daynum, strjoin(Vars)));
end

% The following will be executed for cnmf_e(BatchVer) without the need for
% motion correction
%% 5 "massive" procedure: Extract A from each file
neuron_batch(length(filelist)) = struct('ind_del',[],'signal',[],'FileOrigin',[],'neuron',[]);

parfor i= 1:length(filelist)  
    mode='massive';
    nam=fullfile(datadir,filelist(i).name);    
    %[~,neuron_batch(i)]=demo_endoscope2(gSig,gSiz,min_corr,min_pnr,min_pixel,bd,FS,SSub,TSub,bg_neuron_ratio,nam,mode,[],Afinal,neuron_batch(i),convolveType,merge_thr);
    [~,neuron_batch(i)]=demo_endoscope2(bg_neuron_ratio,merge_thr,with_dendrites,K,sframe,num2read,...
                                   nam,neuron_full,mode,[],neuron_batch(i),Afinal,...
                                   thresh_detecting_frames);
    neuron_batch(i).FileOrigin=filelist(i); % save origin(filelist)
end
fprintf('Massive extraction in each file done.');

%% 6 Save A*C and backsub video

parfor i= 1:length(filelist)
    for j=1:size(neuron_batch(i).neuron.A,2)
        jA=neuron_batch(i).neuron.A(:,j);
        jC=neuron_batch(i).neuron.C(j,:);
        neuron_batch(i).signal(j,:)=median(jA(jA>0)*jC);
    end        
    fprintf('neuron_batch %.0f extraction done\n', i);
end

%% 5.5 deconvolve signal
[~, ~, ~, ~,neuron_batch]=PartTraces(neuron_batch);

%fprintf('First %.0f neurons are successfully deconvolved in each file while those after that are missing in some files\n', sum(~ind_del_final));
fprintf('ALL extractions done.\n');
eval(sprintf('save %sCNMFE_BatchVer.mat %s -v7.3', outputdir, 'neuron_batch'));
fprintf('ALL data saved, check them out!');

%% Making videos
outputdir_video='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/videos/';
d1=File(1).options.d1; d2=File(1).options.d2;
MakingVideos(File,d1,d2,num2str(daynum),outputdir_video)
clear File ACS
MakingVideos([],d1,d2,num2str(daynum),outputdir_video,true,datadir,filelist)
fprintf('ALL videos saved, check them out!');