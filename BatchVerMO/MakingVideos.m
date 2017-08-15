%% Making videos for each day's data
%datafolder='X:\EmilyShijieShared_old\6922_moBatchVer\';
%datafolder='X:\EmilyShijieShared\ProcessedCalciumData\6991FirstFewDaysForBatch\ActuallyUsedInCNMFE\';
type='*20170713*';
d1=300;
d2=400;

video_datalist=dir(fullfile(datafolder,type));
display(length(video_datalist))
Ysignal=[];
for i=1:length(video_datalist)
    Ysignal_tmp=load(fullfile(datafolder,video_datalist(i).name),'Y');
    [Yest, results] = local_background(Ysignal_tmp.Y, 1, 17, [], [], 5);
    clear Ysignal_tmp
    Ysignal=cat(3,Ysignal,Yest);
end

%%
v = VideoWriter([type '.mp4']);
open(v)
bgPermute=Ysignal/1000;
clear Ysignal
img1 = max(bgPermute, 0);
img1 = img1 ./ max(img1(:));
outputVideo=VideoWriter([datafolder type(2:end-1) 'ProcessedSignal.avi']);
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