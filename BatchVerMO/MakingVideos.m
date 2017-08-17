function MakingVideos(File,d1,d2,currentday,outputdir,datadir,filelist,kt)
%% Making videos for each day's data
% kt: avi_file.FrameRate= neuron.Fs/kt;

% some parameters for default
t_begin = 1;

Ysignal=[];
Yac=[];
Ybg=[];
center_ac=[];

for i=1:length(File)
    Ysignal=[Ysignal,File(i).Ysignal];
    Yac=[Yac,File(i).neuron.A*File(i).neuron.C];
    center_ac=[center_ac; max(File(i).neuron.A,[],1)'.*max(File(i).neuron.C,[],2)];
    Ybg=[Ybg,File(i).Ybg];
end
Ysignal=reshape(Ysignal,d1,d2,[]);
Yac=reshape(Yac,d1,d2,[]);
Ybg=reshape(Ybg,d1,d2,[]);
center_ac=median(center_ac);
t_end=size(Yac,3);


Y=[];
for i=1:length(filelist)
    Y_tmp=matfile(fullfile(datadir,filelist(i).name));
    Y=cat(3,Y,Y_tmp.Y);
end
Y=double(Y);

figure('position', [0,0, 600, 400]);


range_res = [-1,1]*center_ac;



if ~exist('range_ac', 'var')
    range_ac = center_ac*1.01+range_res;
end

if ~exist('range_Y', 'var')
    if ~exist('multi_factor', 'var')
        temp = quantile(Y(randi(numel(Y), 10000,1)), [0.01, 0.98]);
        multi_factor = ceil(diff(temp)/diff(range_ac));
    else
        temp = quantile(Y(randi(numel(Y), 10000,1)), 0.01); 
    end
    center_Y = temp(1) + multi_factor*center_ac;
    range_Y = center_Y + range_res*multi_factor;
end
%% create avi file

avi_file = VideoWriter([outputdir currentday '_Videos.avi']);
if ~isnan(File(1).neuron.Fs)
    avi_file.FrameRate= File(1).neuron.Fs/kt;
end
avi_file.open();

ax_y =   axes('position', [0.015, 0.51, 0.3, 0.42]);
ax_bg=   axes('position', [0.015, 0.01, 0.3, 0.42]);
ax_signal=    axes('position', [0.345, 0.51, 0.3, 0.42]);
ax_denoised =    axes('position', [0.345, 0.01, 0.3, 0.42]);
ax_res =    axes('position', [0.675, 0.51, 0.3, 0.42]);

for m=t_begin:kt:t_end
    axes(ax_y); cla;
    imagesc(Ybg(:, :,m)+Ysignal(:, :, m), range_Y);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title('Raw data');
    axis equal off tight;
    
    axes(ax_bg); cla;
    imagesc(Ybg(:, :, m),range_Y);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    axis equal off tight;
    title('Background');
    
    axes(ax_signal); cla;
    imagesc(Ysignal(:, :, m), range_ac); hold on;
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title(sprintf('(Raw-BG) X %d', multi_factor));
    axis equal off tight;
    
    axes(ax_denoised); cla;
    imagesc(Yac(:, :, m), range_ac);
    %     imagesc(Ybg(:, :, m), [-50, 50]);
    title(sprintf('Denoised X %d', multi_factor));
    axis equal off tight;
    
    axes(ax_res); cla;
    imagesc(Ysignal(:, :, m)-Yac(:, :, m), range_res);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title(sprintf('Residual X %d', multi_factor));
    axis equal off tight;
    %         subplot(4,6, [5,6,11,12]+12);
    
    drawnow();
    temp = getframe(gcf);
    temp = imresize(temp.cdata, [400, 600]);
    avi_file.writeVideo(temp);
end

avi_file.close();
end

% bgPermute=Ysignal/1000;
% clear Ysignal
% img1 = max(bgPermute, 0);
% img1 = img1 ./ max(img1(:));
% outputVideo=VideoWriter([outputdir currentday 'ProcessedSignal.avi']);
% outputVideo.FrameRate=30;
% open(outputVideo);
% 
% for p = 1:size(img1,3)
%     img = img1(:,:,p);
%     writeVideo(outputVideo,img);
% end
% 
% close(outputVideo);


%Ysignal=reshape(Ysignal,d1,d2,[]);

% bgPermute=Ysignal/1000;
% clear Ysignal
% img1 = max(bgPermute, 0);
% img1 = img1 ./ max(img1(:));
% outputVideo=VideoWriter([outputdir currentday 'RawSignal.avi']);
% outputVideo.FrameRate=30;
% open(outputVideo);
% 
% for p = 1:size(img1,3)
%     img = img1(:,:,p);
%     writeVideo(outputVideo,img);
% end
% 
% close(outputVideo);
