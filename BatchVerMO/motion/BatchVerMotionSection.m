%% cnmfe (BatchVer) - MotionSection
%  Shijie Gu
%To do: delete days that have no neurons.

%% 1. Concatenate A's from all CNMFE_BatchVer.mat's

outputdir='/Users/gushijie/Documents/Fee/BatchedVerDebug/6991/BatchVerResult/';
%DatadirForA=fullfile(outputdirForA,'*Afinal*');
DatadirForA=fullfile(outputdir,'*Afinal*');
Alist=dir(DatadirForA);
AnumInFolder=numel(Alist);
%%
AsfromDays=[];
AsfromDaysCell={};
sizes=[];
for ia=1:AnumInFolder
    Anext=load(fullfile(outputdir,Alist(ia).name));
    %ColorAllNeurons(Anext.Afinal,300,400,num2str(ia),outputdirForA)
    Atemp=Anext.Afinal;
    AsfromDaysCell{ia}=Atemp;   % individual column is each A, for actual motion correction           
    k = size(Atemp,2);
    sizes=[sizes k];  
    
    C=ones(k,1);    
    Brainbow = Atemp*C; 
    Brainbow = reshape(Brainbow,300,400);  %%%%%%%%%%%%%%%% replace 300,400 with some other
    AsfromDays=cat(3,AsfromDays,Brainbow); % whole picture, for registering, getting shifts        
end

%% 2. for each A in AsfromDays, register others to it.
Y = AsfromDays;
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));

options_nonrigid = NoRMCorreSetParms('upd_template',true,'iter',1,...
                                     'd1',size(Y,1),'d2',size(Y,2),'grid_size',[80,80],'min_patch_size',[50,100,1],'overlap_pre',[10,10,1],'overlap_post',[10,10,1],...
                                     'mot_uf',1,'bin_width',1,...
                                     'max_shift',50,'max_dev',20,'us_fac',5,...
                                     'boundary','zero','iter',1);
                                 % mot_uf: upsamling factor for
                                           % interpolation and individual
                                           % registration.
                                 % us_fac: Upsampling factor (integer). Images will be registered to 
                                           %   within 1/usfac of a pixel. For example usfac = 20 means the
                                           %   images will be registered within 1/20 of a pixel. (default = 1)
%%
[~, ~, ~] = maybe_spawn_workers(4); 
%init_par_rng(2016);
%%
M=cell(1,AnumInFolder);   
shifts=cell(1,AnumInFolder);
for ia=5:5%1:AnumInFolder
    Y_ex_oneday=Y;
    Y_ex_oneday(:,:,ia)=[];
    Y_ex_oneday=Y_ex_oneday(:,:,[4,3,2,1]);
    template_oneday=Y(:,:,ia);
    As_temp=AsfromDaysCell;
    As_temp(ia)=[];
    As_temp=As_temp([4,3,2,1]);
    siz_temp=sizes;
    siz_temp(ia)=[];
    siz_temp=siz_temp([4,3,2,1]);
    tic; [M{ia},shifts{ia},template,xxsfyysf] = normcorre_BatchVer(Y_ex_oneday,options_nonrigid,template_oneday,siz_temp,As_temp,[1,300,61,400,1,1]); toc
    M_temp{1}=M{5}{4};
    M_temp{2}=M{5}{3};
    M_temp{3}=M{5}{2};
    M_temp{4}=M{5}{1};
    M_temp{5}=AsfromDaysCell{ia};
    M{ia}=M_temp;
%     M_temp=M{ia}(ia:end);
%     M{ia}{ia}= AsfromDaysCell{ia};
%     M{ia}(ia+1:AnumInFolder)= M_temp;    
end
%%
for ia=5:5%numel(M)
    for ii=1:numel(M{ia})
        ColorAllNeuronsForMo(M{ia}{ii},300,400,[num2str(ii) '8080UpdatingMondayMor after' num2str(ia)],'/Users/gushijie/Documents/Fee/BatchedVerDebug/6991/BatchVerResult/Cleantry/',xxsfyysf,[10,10]);
    %ColorAllNeuronsForMo(A2,300,400,['18' 'Newupdating after' num2str(ia)],outputdirForA,xxsfyysf);
    end
end




