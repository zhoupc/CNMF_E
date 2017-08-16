function MakingVideos(File,d1,d2,currentday,outputdir,rawOrNot,datadir,filelist)
%% Making videos for each day's data
%datafolder='X:\EmilyShijieShared_old\6922_moBatchVer\';
%datafolder='X:\EmilyShijieShared\ProcessedCalciumData\6991FirstFewDaysForBatch\ActuallyUsedInCNMFE\';


%Version1
%type='*20170713*';
% d1=300;
% d2=400;

% video_datalist=dir(fullfile(datafolder,type));
% display(length(video_datalist))
% Ysignal=[];
% for i=1:length(video_datalist)
%     Ysignal_tmp=load(fullfile(datafolder,video_datalist(i).name),'Y');
%     [Yest, results] = local_background(double(Ysignal_tmp.Y), 1, 17, [], [], 5);
%     clear Ysignal_tmp
%     Ysignal=cat(3,Ysignal,Yest);
% end

%Version2
if nargin<6
    rawOrNot=false;
    filelist=[];
    datadir=[];
end

if rawOrNot==false
    Ysignal_all=[];
    for i=1:length(File)
        Ysignal_all=[Ysignal_all,File(i).Ysignal];
    end
    Ysignal=reshape(Ysignal_all,d1,d2,[]);
    

    bgPermute=Ysignal/1000;
    clear Ysignal
    img1 = max(bgPermute, 0);
    img1 = img1 ./ max(img1(:));
    outputVideo=VideoWriter([outputdir currentday 'ProcessedSignal.avi']);
    outputVideo.FrameRate=30;
    open(outputVideo);
    
    for p = 1:size(img1,3)
        img = img1(:,:,p);
        writeVideo(outputVideo,img);
    end
    
    close(outputVideo);

%
% Ysignal=[];
% for i=1:length(File)
%     Ysignal=cat(4,reshape(File(i).Ysignal,d1,d2,1,size(File(i).Ysignal,2)));
% end
% tmp = bsxfun(@rdivide, Ysignal, max(Ysignal,[],4));
% clear File
% writeVideo(v,tmp)
% close(v)
else
    Ysignal_all=[];
    for i=1:length(filelist)
        Y=matfile(fullfile(datadir,filelist(i).name));
        Ysignal_all=[Ysignal_all,Y.Y];
    end
    Ysignal=reshape(Ysignal_all,d1,d2,[]);
    

    bgPermute=Ysignal/1000;
    clear Ysignal
    img1 = max(bgPermute, 0);
    img1 = img1 ./ max(img1(:));
    outputVideo=VideoWriter([outputdir currentday 'RawSignal.avi']);
    outputVideo.FrameRate=30;
    open(outputVideo);
    
    for p = 1:size(img1,3)
        img = img1(:,:,p);
        writeVideo(outputVideo,img);
    end

    close(outputVideo);
end