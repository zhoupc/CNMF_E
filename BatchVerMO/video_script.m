%% Script Making videos from Ysignal
outputdir='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/';
outputdir_Corr='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/Correlation/';
Filesignal_list=dir(fullfile(outputdir,'*PartI_File*'));
for i=1:length(Filesignal_list)
    nam=Filesignal_list(1).name;
    daynum=nam(1:8);
    load(fullfile(outputdir,['LogisticscnmfeBatchVer' nam(1:8) '.mat']))
    load(fullfile(outputdir,nam))
    
    d1=File(1).options.d1; d2=File(1).options.d2;
    MakingVideos(File,d1,d2,num2str(daynum),outputdir_video)
    clear File
    MakingVideos([],d1,d2,num2str(daynum),outputdir_video,true,datadir,filelist)
    fprintf([daynum 'videos saved, check them out!']);
end
fprintf('ALL videos saved, check them out!');