function CNMFE_Read_In(pathToDataName,OutputFolder)
% CNMFE_Read_In(DataFolder,dataName,outputFolder)
% Example:
% CNMFE_Read_In('/Users/gushijie/Documents/MATLAB/CalciumExtraction/CaELM_row3888.mat','/Users/gushijie/Documents/MATLAB/CNMF_E_New/analysisOut/')
global numFrame num2read ...
     nam saveName outputFolder data nam_mat; %#ok<NUSED> % global variables, don't change them manually
    
nam = pathToDataName;
outputFolder = OutputFolder;


saveName = fullfile(outputFolder, strcat(char(date),' cnmfe_results.mat'));


cur_cd = cd();
if ~exist(outputFolder, 'dir'); 
    mkdir(outputFolder);
else
    fprintf('The folder has been created and old results will be overwritten. \n');
end
cd(outputFolder);
 
cnmfe_choose_data;

globalVars = (who('global'))';
globalVars{1,1+end}='data';
eval(sprintf('save %sCNMFE_Read_In.mat %s -v7.3', outputFolder, strjoin(globalVars)));
fprintf('CNMFEs Read In step Done');
%save(fullfile(outputFolder,'CNMFE_Read_In.mat'), 'data'); 