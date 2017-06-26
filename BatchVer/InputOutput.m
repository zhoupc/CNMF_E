function [datadir,outputdir,kind,running_on_cluster]=InputOutput(varargin)
%  Example:
%  [datadir,outputdir,kind]=InputOutput('datadir','/Volumes/data0-shared/elm/ProcessedCalciumData/7030FirstFewDaysForBatch/JustAFewToTest/','outputdir','/Users/gushijie/Documents/Fee/','datakind','*CaELM*');
%  Input: 1-code directory, 2-data directory, 3-what "kind" of data you want to read in to extract cells.
%         4-running_on_cluster, default is false.
%  "kind" is wildcard in matlab. This is used in dir() function, see
%         documentation dir() for how to use this well.

p = inputParser;
addParameter(p,'datadir',[]);
addParameter(p,'outputdir',[]);
kinddefault='*'; %read in everything in the folder.
addParameter(p,'datakind',kinddefault);
addParameter(p,'running_on_cluster',false); % running on cluster or not
parse(p, varargin{:});

datadir=p.Results.datadir;
outputdir=p.Results.outputdir;
kind=p.Results.datakind;
running_on_cluster=p.Results.running_on_cluster;

Datadir=[datadir,kind];
filelist=dir(Datadir);
numInFolder=numel(filelist);

%first check how many input files are there. Minimum requirement is two.
if numInFolder==0
    ME = MException('cnmfeBatchVer:NotEnoughInput', 'This is an empty folder. \n Please redefine input.');
    throw(ME)
elseif numInFolder==1
    ME2 = MException('cnmfeBatchVer:OnlyOneInput', 'Only one file \n Use normal cnmf-e instead.');
    throw(ME2)
end

% choose final A from every every_file_num as samples in the folders.
prompt=sprintf('%.0f files chosen, \n input x so that every x files to sample as input. Input x then hit "enter": ', numInFolder);
every_file_num = input(prompt);

if isempty(every_file_num)
    every_file_num=max(2,numInFolder/2);
end

choose_ind=mod(1:numel(filelist),every_file_num)==1;

if numel(filelist)<=every_file_num
    error('Zero file for cnmf-e BatchVer sampling. Please lower your every_file_num')
elseif numel(filelist)<2*every_file_num
    error('Only one file for cnmf-e BatchVer sampling. Please lower your every_file_num')
elseif sum(choose_ind)<5
    warning('less than 5 files selected as sample for data. Be aware of this.')
end
    
if running_on_cluster
    [~, ~, ~] = maybe_spawn_workers(4); 
    init_par_rng(2016);
end
end

