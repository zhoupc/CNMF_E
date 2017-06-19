%% Batched CNMF-E version
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

if sum(cellfun('isempty', {filelist.date}))>0
    display(['error in folder:',datadir])
    return
else
    filelist=filelist(mod(1:numel(filelist),3)==1); % pick every other 3 files in the folder as our sample for estimating A.    
end

File(length(filelist)) = struct('HY',[],'A',[],'C',[],'C_raw',[],'S',[],'P',[],'options',[],...
    'Ain',[],'Cin',[],'Cin_raw',[],'STD',[],'Aunique',[],'STDunique',[]);

parfor i= 1:length(filelist)
    nam=fullfile(datadir,filelist(i).name);
    demo_endoscope2
    File(i).A=neuron.A;
    File(i).options=neuron.options;
    File(i).Y=Y
    File(i).Ybg=Ybg; %Ybg+b0
    File(i).Ysignal=Ysignal;
    File(i).Ysignal_sn=Ysignal_sn;
    File(i).noise=noise;    
end
%% 2
A=cat(2,File.A);
Acenter = com(A, d1, d2); %nr x 2 matrix, with the center of mass coordinates
Amask=A>0;
parfor i= 1:length(filelist)
    A2C2A(File, A, Acenter, Amask, i, options)
end
%% 3
merge_thr=[0.7,0.7];
[merged_ROIs,commonA,commonC,ind_del] = mergeAC(Amask,File,merge_thr);


%% 4
parfor i=1:length(filelist)
    FileA=File(i).Ain;
    FileSTD=File(i).STD;
    Aunique(i)=FileA(:,~ind_del);
    STDunique(i)=FileSTD(~ind_del);    
end
weightedA=ReducingA(File);
Afinal=cat(2,commonA,weightedA);

%% 5
filelist=dir(Datadir);
File(length(filelist)) = struct('A',[],'C',[],'C_raw',[],'S',[],'P',[],'options',[]);

parfor i= 1:length(filelist)
    nam=fullfile(datadir,filelist(i).name);
    demo3
    File(i).A=neuron.A;
    File(i).options=neuron.options;
    File(i).Y=Y
    File(i).Ybg=Ybg; %Ybg+b0
    File(i).Ysignal=Ysignal;
    File(i).Ysignal_sn=Ysignal_sn;
    File(i).noise=noise;    
end






    
    
    
    
    
    
    