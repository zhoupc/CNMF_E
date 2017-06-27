%% Demo for cnmfe (BatchVer)
%  Shijie Gu, techel@live.cn

%% A. Input on your PC
%  (1) code directory, data directory, sample directory, what to sample, what to extract.
%  (2) normal CNMF-E parameters.
%  (3) other parameters.

% codeDir='/home/shijiegu/cnmf_e/';
% addpath(genpath(codeDir));
% codeDir2='/home/shijiegu/caprocessing/';
% addpath(genpath(codeDir2));

%(1)
[datadir,sampledir,outputdir,filelist,samplelist]=...
           InputOutput('datadir','/Volumes/data0-shared/elm/ProcessedCalciumData/7030FirstFewDaysForBatch/JustAFewToTest/',...
                       'sampledir',[],...
                       'outputdir','/Users/gushijie/Documents/Fee/',...
                       'datakind','*CaELM*',...
                       'samplekind','*CaELM*',...
                       'SamplingMethod','manual',...
                       'running_on_cluster','true');

%(2)
gSig=10;
gSiz=17;
min_corr=0.85;
min_pnr=6.5;
bg_neuron_ratio = 1;    % spatial range / diameter of neurons
FS=30;
SSub=1;
TSub=1;
namepattern=1:35;       % In sampling stage, for each file running cnmfe, save A's so you can roughly check what neuron is picked in which file.
                        % Each pic's name is from some characters from raw
                        % data's filename. Here I use 1:35
%(3)
% parameters for putting similar neurons in sequence. type help Over_Days_findAnn
correlation_thresh=0.6;
max2max2nd=1.1;
skewnessthresh=0;
% Merge similar neurons based on spatial AND temporal correlation
merge_thr=[0.7,0.7]; 

% Finally, save all these parameters/variables into LogisticscnmfeBatchVer.mat
save([outputdir 'LogisticscnmfeBatchVer.mat'])

%% B. Work on cluster from now
load('LogisticscnmfeBatchVer.mat')
cnmfeBatchVer_ClusterPart