%% Batched CNMF-E version
%% 1
%%%% Input1 %%%%
codeDir='/Users/gushijie/Documents/MATLAB/CaImaging';
addpath(genpath(codeDir));
datadir='/Volumes/data0-shared/elm/ProcessedCalciumData/sleep/7030/061117_5140F/';
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

parfor i= 1:length(filelist)
    picname=filelist(i).name
    nam=fullfile(datadir,filelist(i).name);
    mode='initiation';
    [A0s{i},File(i)]=demo_endoscope2(gSig,gSiz,min_pnr,bg_neuron_ratio,nam,mode,picname,[]); 
end
[ns_storage_1]=Over_Days_findAnn(Ains,0.6,1.1,6);
%% 2
for i= 1:length(filelist)
    A0=A0s{i};
    corr_ind=ns_storage_1(:,1);
    unique_ind=setdiff(1:size(A0,2),corr_ind);
    A0s{i}=[A0(:,corr_ind) A0(:,unique_ind)];
end

A=cat(2,A0s{:});
Acenter = com(A, d1, d2); %nr x 2 matrix, with the center of mass coordinates
Amask=A>0;

ACS(length(filelist)) = struct('Ain',[],'Cin',[],'STD',[]);

parfor i= 1:length(filelist)
    for j=1:length(filelist)
        ACS(i)=A2C2A(File(i), A0s, j, options, 1);
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

%% 5
filelist=dir(Datadir);
%AC(length(filelist)) = struct('Ain',[],'Cin',[]); % AC is used in initiation in demo_2
As=cell(1,length(filelist));
Cs=cell(1,length(filelist));

parfor i= 1:length(filelist)
    
    picname=filelist(i).name
    nam=fullfile(datadir,filelist(i).name);
    mode='massive';
    [As{i},Cs{i}]=demo_endoscope2(gSig,gSiz,min_pnr,bg_neuron_ratio,nam,mode,picname,Afinal);
    
end
%% 6
[ns_storage]=Over_Days_findAnn(As,0.6,1.1,6); % most senstive to skewness border line 5 and 6.

for i=1:3 %length(As)
    subplot(1,length(As),i)
    circleNeuronsFromA(As{i}, ns_storage(:,i))
    title({['neuron from day ' num2str(i)], ['with ' num2str(size(ns_storage,1)) ' neurons tracked.']})
end   