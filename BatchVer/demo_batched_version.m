%% Batched CNMF-E version
%  Shijie Gu, techel@live.cn

%  Direction of use:
%  Input in Section1, let other sections run automatically.
%% 1 Input 
%  (1) code directory, data directory, how many every files to sample.
%  (2) normal CNMF-E parameters.

codeDir='C:\Users\emackev\Documents\MATLAB\cnmf_e';
addpath(genpath(codeDir));
codeDir2='C:\Users\emackev\Documents\MATLAB\CaProcessing';
addpath(genpath(codeDir2));
datadir='\\feevault\data0\ProcessedCalciumData\sleep\7030\061117_5140F\';
kind='*CaELM*';
Datadir=[datadir,kind];
filelist=dir(Datadir);
every_file_num=7;       % choose final A from every every_file_num as samples in the folders.

running_on_cluster=true; % running on cluster or not

gSig=10;
gSiz=17;
min_pnr=6.5;
bg_neuron_ratio = 1;    % spatial range / diameter of neurons
% also "picname" in section2 should be catered to your way of naming files.

if running_on_cluster
    [poolObj, ownPool, poolSize] = maybe_spawn_workers(4); 
    init_par_rng(2016);
end
%% 2 
%%% Preparation for getting A from each sample file.
if sum(cellfun('isempty', {filelist.date}))>0
    display(['error in folder, might be an empty folder:',datadir])
    return
else
    filelist=filelist(mod(1:numel(filelist),every_file_num)==1); % pick every some files in the folder as our sample for estimating A.    
end

% pre-allocate for parfor loop. 
%File is the main output structure of demo_endoscope2. It is different for
%sampling process (mode="initiation") and for running for each file with a given
%A, (mode="massive")
    % containing each file's Y or Ysignal which is denoised and background
    % subtracted data.
%A0s is the output structure containing all A's obtained from each of the
    %sample files.
    
File(length(filelist)) = struct('options',[],'Y',[],'Ysignal',[]);
A0s=cell(1,length(filelist));

%%% Running normal CNMF-E for each file
parfor i= 1:length(filelist)
    picname=filelist(i).name(1:35) % For each file, save A's so you can roughly check what neuron is picked in which file. 
    nam=fullfile(datadir,filelist(i).name);
                                    % The new cnmf-e demo is converted to
                                    % function to cater to parfor loop.
                                    % Meanwhile, some part of cnmf-e demo
                                    % is applied to this sampling process
                                    % while some simplied methods are
                                    % applied for later for each file.
                                    % Use 'mode' variable to specify this.
    mode='initiation';
    [A0s{i},File(i)]=demo_endoscope2(gSig,gSiz,min_pnr,bg_neuron_ratio,nam,mode,picname,[],File(i));
    fprintf('Sampling file number %.0f done\n', i);
end
%%
%%% Order similar neurons in the same sequence in each file, not necessary,
%%% but nice to do. It is fast.
[ns_storage_1]=Over_Days_findAnn(A0s,0.6,1.1,0);
for i= 1:length(filelist)
    A0=A0s{i};
    corr_ind=ns_storage_1(:,1);
    unique_ind=setdiff(1:size(A0,2),corr_ind);
    A0s{i}=[A0(:,corr_ind) A0(:,unique_ind)];
end

%%%  This A is raw and plain, it is just all A's from all files concatnated
%%%  together.
A=cat(2,A0s{:});
Amask=A>0;

%%% Use this A, in each file i, find C's corresponding to each A's found in
%%% file j.
ACS(length(filelist)) = struct('Ain',[],'Cin',[],'STD',[]);

parfor i= 1:length(filelist)
    for j=1:length(filelist)
        ACS(i)=A2C2A(ACS(i),File(i), A0s, j, File(i).options, 1);
    end
end
%% 3 Determine commonA's
%%% Merge similar neurons based on spatial AND temporal correlation
merge_thr=[0.7,0.7];
[merged_ROIs,commonA,commonC,ind_del] = mergeAC(Amask,ACS,merge_thr);

%% 4 Determine uniqueA's
Aunique=cell(1,length(filelist));
STDunique=cell(1,length(filelist));
parfor i=1:length(filelist)
    FileA=ACS(i).Ain;
    FileSTD=ACS(i).STD;
    Aunique{i}=FileA(:,~ind_del);
    STDunique{i}=FileSTD(~ind_del);    
end
weightedA=ReducingA(Aunique,STDunique);

%% 5 Determine Afinal that will be used to extract C's in each file.
Afinal=cat(2,commonA,weightedA);

%%% Some processes making Afinal nicer, modified from Pengcheng Zhou's
%%% idea.
for i=1:size(Afinal,2)
    ai=Afinal(:,i);
    temp = full(ai>quantile(ai, 0.5, 1));
%     l = bwlabel(reshape(temp, File(1).options.d1, File(1).options.d2), 4); 
%     temp(l~=l(ind_ctr)) = false; 
    ai(~temp(:)) = 0;
    Afinal(:,i)=ai;
end

nz_ind=any(Afinal);
Afinal=Afinal(:,nz_ind);

%% 6 "massive" procedure: Extract A from each file
filelist=dir(Datadir);
FILE(length(filelist)) = struct('A',[],'C',[],'ind_del',[]);
parfor i= 1:2 %length(filelist)  
    picname=filelist(i).name(1:35)
    nam=fullfile(datadir,filelist(i).name);
    mode='massive';
    [~,FILE(i)]=demo_endoscope2(gSig,gSiz,min_pnr,bg_neuron_ratio,nam,mode,picname,Afinal,[]);    
end
fprintf('Main extraction done');

neuron(length(filelist)) = struct('signal',[],'filelist',[]);
%%% Partition between those neurons found in each file and those not.
ind_del_final_cat=cat(2,FILE.ind_del);
ind_del_final=any(ind_del_final_cat,2);
parfor i= 1:2 %length(filelist)
    nA=FILE(i).A;
    FILE(i).A=[nA(:,~ind_del_final) nA(:,ind_del_final)];
    fprintf('A extraction done');
    nC=FILE(i).C;
    FILE(i).C=[nC(~ind_del_final,:);nC(ind_del_final,:)];
    fprintf('C extraction done');
    %%% save this data's origin

    for j=1:size(FILE(i).A,2)
        jA=FILE(i).A(:,j);
        jC=FILE(i).C(j,:);
        neuron(i).signal(j,:)=median(jA(jA>0)*jC);
        neuron(i).filelist=filelist(i);
    end        
    fprintf('FILE %.0f extraction done\n', i);
end
fprintf('First %.0f neurons are found in each files while those after that are missing in some files', sum(ind_del_final));
fprintf('ALL extraction done');

%% 7 Post-analysis

