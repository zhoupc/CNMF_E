%% First Part (3h)
codeDir='/home/shijiegu/CNMF_E_New';
display('Dir done')
addpath(genpath(codeDir));
path2data='/home/shijiegu/feevault/ProcessedCalciumData/6938_FirstTutNewSyll/compiled.mat';
path2out='/om/user/shijiegu/analysisOut_April05/';

% no input
CNMFE_Read_In(path2data,path2out)

% input needed
CNMFE_Basic_Parameters(15,25,false,'ar1','thresholded',true,true,true);

global Ysiz
patches = construct_patches(Ysiz(1:2),[64 64],[25 25],[16 16]);
[RESULTS] = CNMFE_Load_PatchY(patches, 1, 15676,1, 1);