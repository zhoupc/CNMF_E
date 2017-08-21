%% Making correlation images from Y
function making_video_batch(cnmfefolder,datadir,d1,d2)
%outputdir='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/';

outputdir_video=[cnmfefolder 'videos/'];

%
if ~exist(outputdir_video,'dir')
    mkdir(outputdir_video)
end
cd(outputdir_video)

Filesignal_list=dir(fullfile(cnmfefolder,'*PartI_File*'));
log_list=dir(fullfile(cnmfefolder,'*LogisticscnmfeBatchVer*'));
for i=1:numel(Filesignal_list) %go through days
    File_temponeday=load(fullfile(cnmfefolder,Filesignal_list(i).name));
    lognam=log_list(i).name;
    fprintf(['Working on file ' lognam]);
    load(fullfile(cnmfefolder,lognam))
    MakingVideos(File_temponeday.File,d1,d2,lognam(end-7:end),outputdir_video,datadir,filelist,1)
end
fprintf('ALL videos saved, check them out!');