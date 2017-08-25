%% Making correlation images from Y
function making_video_batch(cnmfefolder,datadir,d1,d2)
    % General description: Making raw and background subtraced and denoised videos of each day's
    %   file. The actual data concatenation and plotting are in the
    %   function MakingVideos.
    % Input: cnmfefolder: In this folder, there must be files for each day
    %                     that has Ysignal, Ybg, neuronA, neuronC.
    %                     The current version uses data in these files:
    %                     "PartI_File","PartI_Afinalsam" for each day, and
    %                     all days' "neuron_batchMO".
    %        datadir: where Y, raw signal files are.
    %        d1 and d2 for reshaping data.
    
    % Shijie Gu, modified from cnmfe_save_video in CNMF-E by Pengcheng Zhou.

    % make folders for videos if necessary.
    outputdir_video=[cnmfefolder 'videos/'];
    if ~exist(outputdir_video,'dir')
        mkdir(outputdir_video)
    end
    cd(outputdir_video)

    Filesignal_list=dir(fullfile(cnmfefolder,'*PartI_File*'));
    samplelist_list=dir(fullfile(cnmfefolder,'*PartI_Afinalsam*'));
    load(fullfile(cnmfefolder,'neuron_batchMO'))
    
    str_p=0; % current day's data's starting point in the neuron_batchMO structure;
    for i=1:numel(Filesignal_list) %go through days
        File_temponeday=load(fullfile(cnmfefolder,Filesignal_list(i).name));
        samplelist_temponeday=load(fullfile(cnmfefolder,samplelist_list(i).name),'samplelist');
        samplelist=samplelist_temponeday.samplelist;
        lognam=Filesignal_list(i).name;
        fprintf(['Working on day ' lognam]);
        % Input:
        % - data source: File, neuron_batchMO(with reference of str_p),datadir+samplelist,
        % - some parameters: d1, d2.
        % - filename
        % - output directory
        MakingVideos(File_temponeday.File,neuron_batchMO,str_p,datadir,samplelist,d1,d2,lognam(1:8),outputdir_video,1)
        str_p=str_p+length(File_temponeday.File); %update new starting point
    end
fprintf('ALL videos saved, check them out!');