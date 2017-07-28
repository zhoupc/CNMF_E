%% Demo for cnmfe (BatchVerMO)
%  modified by Shijie Gu and Emily Mackevicius

%% A. Input on your PC
%  (1) code directory, data directory, sample directory, what to sample, what to extract.
%  (2) normal CNMF-E parameters.
%  (3) other parameters.
clear all

% I TOTAL DAYS, CODE and OUTPUT DIRECTORY
%codedir=...
%addpath(genpath(codedir));
codeDir='/home/shijiegu/cnmf_e/'; % codeDir2='/home/shijiegu/caprocessing/';

outputdir='/Volumes/shared/EmilyShijieShared_old/6922_moBatchVer/';
outputdir_local=outputdir;
outputdir=strrep(outputdir,'/Volumes/shared/','/net/feevault/data0/shared/');
if ~exist(outputdir_local,'dir')
    mkdir(outputdir_local)
end

% II PARAMETERS
% (1) Normal CNMF-E parameters
gSig=10;
gSiz=17;
min_corr=0.8;
min_pnr=9;
min_pixel = 50;     % minimum number of nonzero pixels for each neuron
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

% (2) parameters for putting similar neurons in sequence.
correlation_thresh=0.6; % Those neuron pairs with cross correlation coefficient above 
                        % correlation_thresh will be considered similar, and
                        % will be put in similar sequence in initiation.
max2max2nd=1.1;         % type help Over_Days_findAnn
skewnessthresh=0;       % type help findn1n2
% Merge similar neurons based on spatial AND temporal correlation
merge_thr_2=[0.7,0.5];

%(4) Cluster paramter: worker number
running_on_cluster=true;
workersnum=4;
                     

% III Specify where data are on your local machine.
totaldays=[20170712,20170713,20170714,20170715,20170716,20170717,20170718,20170719];
for i=1:length(totaldays)    
    daynum=totaldays(i);    
    outputdirDetails = [outputdir, num2str(totaldays(i)), '/'];
    
    [datadir,sampledir,~,filelist,samplelist]=...
           InputOutput('datadir','/Volumes/shared/EmilyShijieShared/ProcessedCalciumData/6991FirstFewDaysForBatch/',...
                       'datakind',['*' num2str(totaldays(i)) '*']);
    % replace for actual dir on cluster
    datadir=strrep(datadir,'/Volumes/shared/','/net/feevault/data0/shared/');
    sampledir=strrep(sampledir,'/Volumes/shared/','/net/feevault/data0/shared/');
    
   % Finally, save all these parameters/variables into LogisticscnmfeBatchVer.mat
    ininame='LogisticscnmfeBatchVer';
    save([outputdir_local ininame num2str(totaldays(i)) '.mat'])     
end

%% B. making cluster script
jobdir='/home/shijiegu/jobs/6922BatchVerMo/';
local_jobdir='/Users/gushijie/Documents/MATLAB/Jobs/6922BatchVerMo/';

if ~exist(local_jobdir,'dir')
    mkdir(local_jobdir)
end

cd(local_jobdir)

shellname='BatchVerMO6922.sh';

for i=1:length(totaldays)
    scriptname=sprintf('BatchVerMO6922_%.0f.sbatch',totaldays(i));
    fileID = fopen(scriptname,'w');
    request=['#!/bin/bash\n',...
        '\n',...
        '#SBATCH -n 1\n',...     
        '#SBATCH --cpus-per-task=8\n',...
        '#SBATCH --mem=200000\n',...
        '#SBATCH -t 0-10:00\n',...
        '#SBATCH --time-min=0-01:00\n'];
    fprintf(fileID, request);
    
    errorANDoutAndcd= ['#SBATCH -o %sjob_%%A_%s.out\n',...
        '#SBATCH -e %sjob_%%A_%s.err\n',...
        '\n'];
    fprintf(fileID,errorANDoutAndcd,outputdir,num2str(totaldays(i)),outputdir,num2str(totaldays(i)));
    
    MATLABcommands=['module add mit/matlab/2016b\n',...
        'matlab -nodisplay -singleCompThread -r "addpath(genpath(''%s''));\\\n',...
        'load(''%s''); cnmfeMoBatchVer_ClusterPart"'];
    fprintf(fileID,MATLABcommands,codeDir,[outputdir,'LogisticscnmfeBatchVer',num2str(totaldays(i)),'.mat']);
    
    fclose(fileID);
end
%chmod g+rwx /net/feevault/data0/shared/EmilyShijieShared

%% C. making the shell script

% fileID = fopen(shellname,'w');
% shtext=['#!/bin/sh\n',...
%     'cd ''%s''\n',...
%     'LogisticsLIST=($(ls -1 *Logistics*))\n',...
%     'declare -p LogisticsLIST > .bash_arrays\n',...
%     'NUMLogistics=${#LogisticsLIST[@]}\n',...
%     'ZBNUMLogistics=$(($NUMLogistics - 1))\n',...
%     'echo $ZBNUMLogistics\n',...
%     'sbatch --array=0-$ZBNUMLogistics ''%s'''];
% fprintf(fileID,shtext,outputdir,fullfile(jobdir,scriptname));
% fclose(fileID);
