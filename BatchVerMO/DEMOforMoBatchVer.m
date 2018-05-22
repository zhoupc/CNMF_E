%% Demo for cnmfe (BatchVerMO)
%  modified by Shijie Gu and Emily Mackevicius

%% A. Input on your PC
%  I normal CNMF-E parameters and other new parameters for CNMF-E.
%  II 1/code directory...  [codeDir]
%     2/output directory...[outputdir,outputdirDetails,datadir,filelist,samplelist]
%     3/data directory...sample directory, what to sample, what to extract.
%     4/days               [totaldays]

clear all
% I PARAMETERS
% ----(2)
% about loading data
sframe=1;                                          % first frame to read (optional, default:1)
num2read=[];                                       % how many frames to read   (optional, default: until the end)

% create Source2D class object for storing results and parameters  
with_dendrites=false;                              % with dendrites or not
K = [];                                            % maximum number of neurons to search. 
                                                   % Use [] to search the number automatically
neuron_full = Sources2D('ssub',1, 'tsub',1, ...   % downsampling
    'gSig',10,...                                 % sigma of the 2D gaussian (width of the gaussian kernel) that approximates cell bodies
    'gSiz',17,...                                 % maximum diameter of neurons in the image plane. larger values are preferred.
    'min_corr',0.8,'min_pnr',14, ...
    'min_pixel',50,...                            % minimum number of nonzero pixels for each neuron
    'bd',1,...                                    % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
    'center_psf',true);                          
neuron_full.Fs = 30;                               % frame rate  
% options for running deconvolution 
neuron_full.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    ...% convolveType: string, defines the model of the deconvolution kernel. possible options are:
    ...                                     % 'ar1':  auto-regressive model with order p=1
    ...                                     % 'ar2':  auto-regressive model with order p=2
    ...                                     % 'exp2': the convolution kernel is modeled as the difference of two
    ...                                     %         exponential functions -
    ...                                     %               h(t) = (exp(-t/tau_d) - exp(-t/tau_r)) / (tau_d-tau_r)
    ...                                     % 'kernel':   a vector of the convolution kernel
    'method', 'constrained', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...   % optimize AR coefficients
    'optimize_b', true, ...      % optimize the baseline
    'optimize_smin', true);      % optimize the threshold

% Parameters in updating
bg_neuron_ratio = 1;    % spatial range / diameter of neurons

% parameters for putting similar neurons in sequence. 
correlation_thresh=0.6; % Those neuron pairs with cross correlation coefficient above 
                        % correlation_thresh will be considered similar, and
                        % will be put in similar sequence in initiation.
max2max2nd=1.1;         % type help Over_Days_findAnn
skewnessthresh=0;       % type help findn1n2

% Merge similar neurons based on spatial AND temporal correlation
% "merge_thr_2" here will be used in MergeAC, the BatchVer specific merging, 
% while "merge_thr" in section2 is used for for normal CNMF-E neuron merging in initiation and iteration process.
merge_thr=[0.7,0.5,0];
merge_thr_2=[0.7,0.5];

% Merge neurons that are very close.
dmin=5;

thresh_detecting_frames=10;  % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

%%%% There are a few more techinical parametrs in demo_endoscope2.

% Cluster-specific parameters, depend on the cluster you use and the CPU number you request. 
    % With this workersnum, cnmfeBatchVer_ClusterPart has functions written to create a cluster and a parpool within it. 
    % The design caters to parpool in cluster use. The functions for these are in a seperate folder "HandlingParforbyGalen" which is written by Galen Lynch.
running_on_cluster=true;
workersnum=4;

% II code directory and output directory

Version='BatchVer';               %with motion: MoBatchVer; without motion: BatchVer.
% ---- /1/----------------------------------------------------------------
codeDir='/home/shijiegu/cnmf_e/';   %where code is on the cluster
% ------------------------------------------------------------------------

% ---- /2/Output directory, local and on cluster--------------------------
outputdir_local='/Volumes/shared/EmilyShijieShared/CaExtraction/6938JustSinging/';
outputdir='/net/feevault/data0/shared/EmilyShijieShared/CaExtraction/6938JustSinging/';
%outputdir=strrep(outputdir_local,'/Volumes/data0-shared/elm/ProcessedCalciumData/','/net/feevault/data0/shared/EmilyShijieShared/CaExtraction/6701/');
% ------------------------------------------------------------------------

% -----------------------------ignore this auto area----------------------
if ~exist(outputdir_local,'dir')
    mkdir(outputdir_local)
end
% ------------------------------------------------------------------------

% ---- /3/ "Motion batch", within batch, no motion, between batch there can be motion.-
% If using BatchVer, you can just input one element in totaldays.
%totaldays=6938;
totaldays_act=[6938];
totaldays=6938; %[0612 0613 0614 0615 0617 0618 0620 0621];
% -------------------------------------------------------------------------

% -----/4/ Fill in datadir on local machine, datakind----------------------
for i=1:length(totaldays_act)    
    daynum=totaldays_act(i);
    if strcmp(Version,'MoBatchVer')
        outputdirDetails = [outputdir, num2str(totaldays_act(i)), '/'];
    elseif strcmp(Version,'BatchVer')
        outputdirDetails = [outputdir, 'details', '/'];
    end
    
    [datadir,sampledir,~,filelist,samplelist]=...
           InputOutput('datadir','/Volumes/shared/EmilyShijieShared/ProcessedCalciumData/6938FirstFewDaysForBatch/Singing/',...
                       'datakind',['*' num2str(totaldays_act(i)) '*'],...
                       'DataMethod','auto');
    % if there is any motion within a day
                       %Add another day in totaldays, use ('DataMethod','manual');
    % replace for actual dir on cluster
    datadir=strrep(datadir,'/Volumes/','/net/feevault/data0/');
    sampledir=strrep(sampledir,'/Volumes/','/net/feevault/data0/');
    
   % Finally, save all these parameters/variables into LogisticscnmfeBatchVer.mat
    ininame='LogisticscnmfeBatchVer';
    save([outputdir_local ininame num2str(totaldays_act(i)) '.mat'])     
end
% ----------------------------------------------------------------------------------------

%% B. making cluster script
% To better control each day's cnmfe behavior, each day has its own script,
% thus constinuting one job-one day. For all these scripts, use MATLAB to
% write in jobdir folder with the common prefix as variable "shellname".

% (1) jobdir folder, shellname
% ---- /1/ ----------------------------------------------------------------
jobdir='BatchVer6938';
local_jobdir=['/Users/gushijie/Documents/MATLAB/Jobs/' jobdir '/'];
shellname=jobdir;
% -------------------------------------------------------------------------

% -----------------------------ignore this auto area-----------------------
if ~exist(local_jobdir,'dir')
    mkdir(local_jobdir)
end
cd(local_jobdir)
% -------------------------------------------------------------------------

for i=1:length(totaldays)
    scriptname=sprintf([shellname '_%.0f.sbatch'],totaldays(i));
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

%% B2. making the shell script
fileID = fopen([shellname '.sh'],'w');
forloop=['cd %s\n',...
         'for file in %s; do\n',...
         'sbatch %s_${file}.sbatch\n',...
         'sleep 1\n',...
         'done'];
fprintf(fileID,forloop,jobdir,num2str(totaldays),shellname);

%% B3. (MoBatchVer only) for cnmfeMoBatchVer_ClusterPart2
scriptname=sprintf([shellname '_part2.sbatch']);
fileID = fopen(scriptname,'w');
request=['#!/bin/bash\n',...
    '\n',...
    '#SBATCH -n 1\n',...
    '#SBATCH --cpus-per-task=8\n',...
    '#SBATCH --mem=200000\n',...
    '#SBATCH -t 0-20:00\n',...
    '#SBATCH --time-min=0-01:00\n'];
fprintf(fileID, request);

errorANDoutAndcd= ['#SBATCH -o %sjob_%%A.out\n',...
    '#SBATCH -e %sjob_%%A.err\n',...
    '\n'];
fprintf(fileID,errorANDoutAndcd,outputdir,outputdir);

MATLABcommands=['cd %s\n',...
    'module add mit/matlab/2016b\n',...
    'matlab -nodisplay -singleCompThread -r "addpath(genpath(''%s''));addpath(genpath(''/home/shijiegu/caprocessing/''));\\\n',...
    'outputdir=''%s'';',...
    'load(fullfile(outputdir,''%s''));load(fullfile(outputdir,''cnmfe_BatchVer_PartII_MotionCorrection.mat''));\\\n',...
    'cnmfeMoBatchVer_ClusterPart2"'];
fprintf(fileID,MATLABcommands,outputdir,codeDir,outputdir,['LogisticscnmfeBatchVer',num2str(totaldays(1)),'.mat']);

fclose(fileID);
