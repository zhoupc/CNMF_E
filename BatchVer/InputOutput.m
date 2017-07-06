function [datadir,sampledir,outputdir,filelist,samplelist]=InputOutput(varargin)
%  [datadir,sampledir,outputdir,filelist,samplelist]=InputOutput(varargin)
%  Example: [datadir,sampledir,outputdir,filelist,samplelist]=
%           InputOutput('datadir','/Volumes/data0-shared/elm/ProcessedCalciumData/7030FirstFewDaysForBatch/JustAFewToTest/',
%                       'sampledir','[]'
%                       'outputdir','/Users/gushijie/Documents/Fee/',
%                       'datakind','*CaELM*',
%                       'samplekind','*CaELM*',
%                       'SamplingMethod','manual'
%                       'DataMethod','manual');
%  Input: 1-data directory, 2-sample directory, which can be empty. If
%         'SamplingMethod' is 'auto', the sample directory will be the same as data directory
%         3-output directory, 4/5-what "kind" of data you want to read in to extract cells (for data and sample).
%         6/7-how to choose data for sampling or extraction, either "auto", every some files will be asked for you
%         as an input to sample data; or, "manual", a window will pop up for you to choose data from the sampledir you input.
%  More notes:
%  "kind" is wildcard in matlab. This is used in dir() function, see
%         documentation dir() for how to use this well.
%  When choosing files "manual"ly, select multiple files by 
%         holding down the Shift or Ctrl key and clicking on additional file names.
%         IMPORTANT: All samples must come from one folder!

%  Shijie Gu, techel@live.cn

p = inputParser;
addParameter(p,'datadir',[]);
addParameter(p,'sampledir',[]);
addParameter(p,'outputdir',[]);
kinddefault='*'; %read in everything in the folder.
addParameter(p,'datakind',kinddefault);
addParameter(p,'samplekind',kinddefault);
expectedmethod = {'auto','manual'};
addParameter(p,'SamplingMethod','auto', @(x) any(validatestring(x,expectedmethod)));
addParameter(p,'DataMethod','manual', @(x) any(validatestring(x,expectedmethod)));

parse(p, varargin{:});

datadir=p.Results.datadir;
sampledir=p.Results.sampledir;
outputdir=p.Results.outputdir;
datakind=p.Results.datakind;
samplekind=p.Results.samplekind;
SamplingMethod=p.Results.SamplingMethod;
DataMethod=p.Results.DataMethod;

if strcmp(DataMethod,'auto')
    Datadir=fullfile(datadir,datakind);
    filelist=dir(Datadir);
    numInFolder=numel(filelist);
    %first check how many input data files are there. Minimum requirement is two.
    if numInFolder==0
        ME = MException('cnmfeBatchVer:NotEnoughInput', 'This is an empty data folder or a incorrectly specified folder. \n Please redefine input.');
        throw(ME)
    elseif numInFolder==1
        ME2 = MException('cnmfeBatchVer:OnlyOneInput', 'Only one file \n Use normal cnmf-e instead.');
        throw(ME2)
    end
else
    [FileName,PathName,~] = uigetfile({'*.*';datakind},'cnmfe(BatchVer)-Manually Choose Data',datadir,'MultiSelect','on');
    datadir=PathName; %Just in case in practice, user changed directory in the folder.
    numInFolder=numel(FileName);
        %first check how many input data files are there. Minimum requirement is two.
    if numInFolder==0
        ME = MException('cnmfeBatchVer:NotEnoughInput', 'This is an empty data folder or a incorrectly specified folder. \n Please redefine input.');
        throw(ME)
    elseif numInFolder==1
        ME2 = MException('cnmfeBatchVer:OnlyOneInput', 'Only one file \n Use normal cnmf-e instead.');
        throw(ME2)
    end
    sprintf('Chosen %.0f files to extract in the end. \n', numInFolder);
    for i=1:length(FileName)
        filelist(i)= dir(fullfile(PathName,FileName{i}));
    end
end

if strcmp(SamplingMethod,'auto')
    sampleInFolder=numel(filelist);
    if sampleInFolder==0
        ME3 = MException('cnmfeBatchVer:NotEnoughSampleInput', 'There is no sample input. \n Please redefine sample input.');
        throw(ME3)
    end    
    % choose final A from every every_file_num as samples in the folders.
    prompt=sprintf('%.0f files in the sample folder, \n input x so that every x files to sample as input. Input x then hit "enter": ', sampleInFolder);
    every_file_num = input(prompt);
    if isempty(every_file_num)
        every_file_num=max(2,numInFolder/2);
    end
    choose_ind=mod(1:numel(filelist),every_file_num)==1;
    if numel(filelist)<=every_file_num
        error('Zero file for cnmf-e BatchVer sampling. Please lower your every_file_num')
    elseif numel(filelist)<2*every_file_num
        warning('Only one file for cnmf-e BatchVer sampling. Lowering your every_file_num is suggested.')
    elseif sum(choose_ind)<5
        warning('less than 5 files selected as sample for data. Be aware of this.')
    end
    samplelist=filelist(choose_ind);
else
    if isempty(sampledir)
        sampledir=datadir;
    end
    [FileName,PathName,~] = uigetfile({'*.*';samplekind},'cnmfe(BatchVer)-Manually Choose Sample',sampledir,'MultiSelect','on');
    sprintf('Chosen %.0f files as samples. \n', length(FileName));
    if numel(FileName)==1
        warning('Only one file for cnmf-e BatchVer sampling. Lowering your every_file_num is suggested.')
    elseif numel(FileName)<5
        warning('less than 5 files selected as sample for data. Be aware of this.')        
    end
    sampledir=PathName; %Just in case in practice, user changed directory in the folder.
    samplelist = struct('name',FileName);
end

end

