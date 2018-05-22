%% Demo for cnmfe (BatchVer)
%  modified by Shijie Gu and Emily Mackevicius

%% A. Input on your PC
%  (1) code directory, data directory, sample directory, what to sample, what to extract.
%  (2) normal CNMF-E parameters.
%  (3) other parameters.
clear all
%(0)
%codedir=...
%addpath(genpath(codedir));
codeDir='/home/shijiegu/cnmf_e/'; % codeDir2='/home/shijiegu/caprocessing/';

%(1) Specify where data are on your local machine. This is a parsing
%step.

[datadir,sampledir,outputdir,filelist,samplelist]=...
           InputOutput('datadir','/Volumes/shared/EmilyShijieShared/ProcessedCalciumData/6991FirstFewDaysForBatch/',...
                       'sampledir',[],...
                       'outputdir','/Volumes/shared/EmilyShijieShared/ProcessedCalciumData/6991FirstFewDaysForBatch/Singing-BatchVerResult-071607170718/',...
                       'datakind','*Ca*',...
                       'samplekind','*',...
                       'SamplingMethod','auto',...
                       'dataMethod','manual');
%(1.1, replace for actual dir on cluster)
datadir=strrep(datadir,'/Volumes/shared/','/net/feevault/data0/shared/');
sampledir=strrep(sampledir,'/Volumes/shared/','/net/feevault/data0/shared/');
outputdir_local=outputdir;
outputdir=strrep(outputdir,'/Volumes/shared/','/net/feevault/data0/shared/');
if ~exist(outputdir_local,'dir')
    mkdir(outputdir_local)
end
daynum=201704270501;
Aoutputdir=outputdir;

%(2)
gSig=10;
gSiz=17;
min_corr=0.75;
min_pnr=6;
min_pixel = 50;      % minimum number of nonzero pixels for each neuron
bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)

bg_neuron_ratio = 1;    % spatial range / diameter of neurons
FS=30;
SSub=1;
TSub=1;
convolveType='ar1';    % convolveType: string, defines the model of the deconvolution kernel. possible options are:
                                        % 'ar1':  auto-regressive model with order p=1
                                        % 'ar2':  auto-regressive model with order p=2
                                        % 'exp2': the convolution kernel is modeled as the difference of two
                                        %         exponential functions -
                                        %               h(t) = (exp(-t/tau_d) - exp(-t/tau_r)) / (tau_d-tau_r)
                                        % 'kernel':   a vector of the convolution kernel
merge_thr=[0.7,0.5,0];
%namepattern=:;       % In sampling stage, for each file running cnmfe, save A's so you can roughly check what neuron is picked in which file.
                        % Each pic's name is from some characters from raw
                        % data's filename. Here I use 1:35
%(3)
% parameters for putting similar neurons in sequence. 
correlation_thresh=0.6; % Those neuron pairs with cross correlation coefficient above 
                        % correlation_thresh will be considered similar, and
                        % will be put in similar sequence in initiation.
max2max2nd=1.1;         % type help Over_Days_findAnn
skewnessthresh=0;       % type help findn1n2
% Merge similar neurons based on spatial AND temporal correlation
merge_thr_2=[0.7,0.5];

%(4)
running_on_cluster=true;
workersnum=4;

% Finally, save all these parameters/variables into LogisticscnmfeBatchVer.mat
ininame='LogisticscnmfeBatchVer.mat';
save([outputdir_local ininame])

%% B. making cluster script
fileID = fopen('BatchVerSLURM6991Sing161718.sh','w');
request={'#!/bin/bash'
    ''
    '#SBATCH -n 1'
    '#SBATCH --cpus-per-task=8'
    '#SBATCH --mem=200000'
    '#SBATCH -t 0-10:00'
    '#SBATCH --time-min=0-01:00'};
fprintf(fileID,'%s\n',request{:});

errorANDoutAndcd= ['#SBATCH -o %sjob_%%A.out\n#SBATCH -e %sjob_%%A.err\n',...
    'cd %s\n'];
fprintf(fileID,errorANDoutAndcd,outputdir,outputdir,outputdir);

MATLABcommands=['module add mit/matlab/2016b\n',...
    'matlab -nodisplay -singleCompThread -r "addpath(genpath(''%s''));\\\n',...
    'load(fullfile(''%s'',''%s'')); cnmfeBatchVer_ClusterPart"'];
fprintf(fileID,MATLABcommands,codeDir,outputdir,ininame);

fclose(fileID);
%chmod g+rwx /net/feevault/data0/shared/EmilyShijieShared
