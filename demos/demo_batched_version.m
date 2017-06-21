%% Batched CNMF-E version
%% 1
%%%% Input1 %%%%
codeDir='C:\Users\emackev\Documents\MATLAB\cnmf_e';
addpath(genpath(codeDir));
codeDir2='C:\Users\emackev\Documents\MATLAB\CaProcessing';
addpath(genpath(codeDir2));

datadir='U:\ProcessedCalciumData\sleep\7030\061117_5140F\';
kind='*CaELM*';
gSig=10;
gSiz=17;
min_pnr=6.5;
bg_neuron_ratio = 1;  % spatial range / diameter of neurons

Datadir=[datadir,kind];
filelist=dir(Datadir);

if sum(cellfun('isempty', {filelist.date}))>0
    display(['error in folder:',datadir])
    return
else
    filelist=filelist(mod(1:numel(filelist),7)==1); % pick every other 3 files in the folder as our sample for estimating A.    
end

File(length(filelist)) = struct('options',[],'Y',[],'Ybg',[],'Ysignal',[],...
    'Ysignal_sn',[],'noise',[],...
    'Ain',[],'Cin',[],'STD',[]);
A0s=cell(1,length(filelist));
%%
parfor i= 1:length(filelist)
    picname=filelist(i).name(1:35)
    nam=fullfile(datadir,filelist(i).name);
    mode='initiation';
    [A0s{i},File(i)]=demo_endoscope2(gSig,gSiz,min_pnr,bg_neuron_ratio,nam,mode,picname,[],File(i));
    fprintf('File %.0f done\n', i);
end
[ns_storage_1]=Over_Days_findAnn(A0s,0.6,1.1,0);
%% 2
for i= 1:length(filelist)
    A0=A0s{i};
    corr_ind=ns_storage_1(:,1);
    unique_ind=setdiff(1:size(A0,2),corr_ind);
    A0s{i}=[A0(:,corr_ind) A0(:,unique_ind)];
end
%%
A=cat(2,A0s{:});
Amask=A>0;

ACS(length(filelist)) = struct('Ain',[],'Cin',[],'STD',[]);

parfor i= 1:length(filelist)
    for j=1:length(filelist)
        ACS(i)=A2C2A(File(i), A0s, j, File(i).options, 1);
    end
end
%% 3
merge_thr=[0.7,0.7];
[merged_ROIs,commonA,commonC,ind_del] = mergeAC(Amask,File,merge_thr);

%% 4
Aunique=cell(1,length(filelist));
STDunique=cell(1,length(filelist));
parfor i=1:length(filelist)
    FileA=File(i).Ain;
    FileSTD=File(i).STD;
    Aunique{i}=FileA(:,~ind_del);
    STDunique{i}=FileSTD(~ind_del);    
end
weightedA=ReducingA(Aunique,STDunique);
Afinal=cat(2,commonA,weightedA);

for i=1:size(Afinal,2)
    ai=Afinal(:,i);
    temp = full(ai>quantile(ai, 0.5, 1));
    l = bwlabel(reshape(temp, nr, nc), 4); 
    temp(l~=l(ind_ctr)) = false; 
    ai(~temp(:)) = 0;
    Afinal(:,i)=ai;
end

%% 5
filelist=dir(Datadir);
As=cell(1,length(filelist));
Cs=cell(1,length(filelist));

parfor i= 1:length(filelist)    
    picname=filelist(i).name(1:35)
    nam=fullfile(datadir,filelist(i).name);
    mode='massive';
    [~,neuron{i}]=demo_endoscope2(gSig,gSiz,min_pnr,bg_neuron_ratio,nam,mode,picname,Afinal,[]);        
end
%% 6
