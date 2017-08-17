%% Making correlation images from Y
function corr_ForBatch(outputdir,datadir,d1,d2)
%outputdir='/net/feevault/data0/shared/EmilyShijieShared_old/6922_moBatchVerNYVersion/';
outputdir_Corr=[outputdir 'Correlation/'];

if ~exist(outputdir_Corr,'dir')
    mkdir(outputdir_Corr)
end
cd(outputdir_Corr)

log_list=dir(fullfile(outputdir,'*LogisticscnmfeBatchVer*'));
display(length(log_list))
for i=1:length(log_list)
    log_nam=log_list(1).name;
    load(fullfile(outputdir,log_nam))
    
    Y=[];
    for j=1:length(filelist)
        Y=matfile(fullfile(datadir,filelist(j).name));
        Y=cat(3,Y,Y.Y);
    end
    display(size(Y))
    Y=double(reshape(Y,d1*d2,[]));
    % divide data into multiple patches
    patch_sz = [3, 3];
    r0_patch = round(linspace(1, d1, 1+patch_sz(1)));
    c0_patch = round(linspace(1, d2, 1+patch_sz(2)));
    nr_patch = length(r0_patch)-1; 
    nc_patch = length(c0_patch)-1; 
    Cn = zeros(d1, d2);
    PNR = zeros(d1,d2);
    % compute correlation_image patch by patch
    bd = 20;
    sig=5; % thresholding noise by sig*std()
    for mr = 1:nr_patch
        r0 = max(1, r0_patch(mr)-bd); % minimum row index of the patch
        r1 = min(d1, r0_patch(mr+1)+bd-1); % maximum row index of the patch
        for mc = 1:nc_patch
            c0 = max(1, c0_patch(mc)-bd); % minimum column index of the patch
            c1 = min(d2, c0_patch(mc+1)+bd-1); % maximum column index of the patch
        
            % take the patch from the raw data
            nrows = (r1-r0+1);  % number of rows in the patch
            ncols = (c1-c0+1);  %number of columns in the patch
            Ypatch = double(Y(r0:r1, c0:c1, :));
        
            % spatially filter the data
            %HY = imfilter(Ypatch, psf, 'replicate');
        
            % copute signal to noise ratio
            %HY = reshape(HY, [], T);
            %HY = bsxfun(@minus, HY, median(HY, 2));
            HY_max=max(Ypatch,[],2);%HY_max = max(HY, [], 2);
            Ysig = get_noise_fft(Ypatch);
%             tmp_PNR = reshape(HY_max./Ysig, nrows, ncols);
%             PNR(r0:r1, c0:c1) = max(PNR(r0:r1, c0:c1), tmp_PNR);
        
        
            %  compute loal correlation
            Ypatch(bsxfun(@lt, Ypatch, Ysig*sig)) = 0;
            tmp_Cn = correlation_image(Ypatch, [1,2], nrows, ncols);
            Cn(r0:r1, c0:c1) = max(Cn(r0:r1, c0:c1), tmp_Cn);   
        end
    end
    figure
    imagesc(Cn, [0, 1]); colorbar;
    axis equal off tight;
    fignam=['correlation image of ' num2str(daynum)];
    title(fignam);
    saveas(gcf,fignam);
end
% Filesignal_list=dir(fullfile(outputdir,'*PartI_File*'));
% for i=1:length(Filesignal_list)
%     nam=Filesignal_list(1).name;
%     daynum=nam(1:8);
%     load(fullfile(outputdir,['LogisticscnmfeBatchVer' nam(1:8) '.mat']))
%     load(fullfile(outputdir,nam))
%     
%     d1=File(1).options.d1; d2=File(1).options.d2;
%     MakingVideos(File,d1,d2,num2str(daynum),outputdir_video)
%     clear File
%     MakingVideos([],d1,d2,num2str(daynum),outputdir_video,true,datadir,filelist)
%     fprintf([daynum 'videos saved, check them out!']);
% end
fprintf('ALL correlation images saved, check them out!');