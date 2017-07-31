%% cnmfe (BatchVer) - MotionSection
%  Shijie Gu
%To do: delete days that have no neurons.

%% 1. Concatenate A's from all CNMFE_BatchVer.mat's
% AsfromDaysPic has pictures
% AsfromDaysCell has A's
birdnum='6922';
cnmfedir='X:\EmilyShijieShared_old\6922_moBatchVer\';
DatadirForA=fullfile(cnmfedir,'*PartI_Afinal*');
Alist=dir(DatadirForA);
AnumInFolder=numel(Alist);

AsfromDaysPic=[];
AsfromDaysCell={};
sizes=[];
for ia=1:AnumInFolder
    Anext=load(fullfile(cnmfedir,Alist(ia).name));
    %ColorAllNeurons(Anext.Afinal,300,400,num2str(ia),outputdirForA)
    Atemp=Anext.Afinal;
    AsfromDaysCell{ia}=Atemp;   % individual column is each A, for actual motion correction           
    k = size(Atemp,2);
    sizes=[sizes k];  
    
    C=ones(k,1);    
    Brainbow = Atemp*C; 
    Brainbow = reshape(Brainbow,300,400);  %%%%%%%%%%%%%%%% replace 300,400 with some other
    AsfromDaysPic=cat(3,AsfromDaysPic,Brainbow); % whole picture, for registering, getting shifts        
end

%% 2. for each A in AsfromDays, register others to it.
Y = AsfromDaysPic;
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));

options_nonrigid = NoRMCorreSetParms('upd_template',true,'iter',1,...
                                     'd1',size(Y,1),'d2',size(Y,2),'grid_size',[80,80],'min_patch_size',[50,50,1],'overlap_pre',[10,10,1],'overlap_post',[10,10,1],...
                                     'mot_uf',10,'bin_width',1,...
                                     'max_shift',50,'max_dev',20,'us_fac',5,...
                                     'boundary','zero','iter',1);
                                 % mot_uf: upsamling factor for
                                           % interpolation and individual
                                           % registration.
                                 % us_fac: Upsampling factor (integer). Images will be registered to 
                                           %   within 1/usfac of a pixel. For example usfac = 20 means the
                                           %   images will be registered within 1/20 of a pixel. (default = 1)
% Additional parameters
update_num=2;

%[~, ~, ~] = maybe_spawn_workers(4); 

%%
global shifts
M=cell(1,AnumInFolder);   
shifts=cell(1,AnumInFolder);
template=cell(1,AnumInFolder);
ind_del=cell(1,AnumInFolder);
for ia=1:AnumInFolder
    M_temp=cell(1,AnumInFolder);
    ind_del_temp=cell(1,AnumInFolder);
    [beforeseq_ind,reseq_ind]=sort(abs([1:AnumInFolder]-ia)); % Closer days are aligned first and then others. Calculate the index here.
    
    Y_ex_oneday=Y(:,:,reseq_ind(2:end));
    Y_oneday=Y(:,:,ia);
    As_ex_oneday=AsfromDaysCell(reseq_ind(2:end));
    siz_ex_oneday=sizes(reseq_ind(2:end));
    
    startendgrid=[1,300,61,400,1,1];
    tic; [M{ia},shifts{ia},template{ia},xxsfyysf,ind_del{ia}] = normcorre_BatchVer(Y_ex_oneday,options_nonrigid,Y_oneday,siz_ex_oneday,As_ex_oneday,startendgrid,update_num); toc
    M_temp(reseq_ind(2:end))=M{ia}(1:end);              M_temp{ia}=AsfromDaysCell{ia};       M{ia}=cell2mat(M_temp); 
    ind_del_temp(reseq_ind(2:end))=ind_del{ia}(1:end);  ind_del_temp{ia}=false(1,sizes(ia)); ind_del{ia}=cell2mat(ind_del_temp);
end
            
%%
ind_del_full=sum(reshape(cell2mat(ind_del),1,sum(sizes),[]),3)>0;    %sum(cat(3,ind_del{:}))>0;
ind_del_full_cell=mat2cell(~ind_del_full,1,sizes);
N_eachday=cellfun(@(x) sum(x), ind_del_full_cell);
M_del = cellfun(@(x) x(:,~ind_del_full), M, 'UniformOutput',0);
M_final = cellfun(@(x) mat2cell(x, [size(x,1)], N_eachday), M_del, 'UniformOutput',0);



%%
for ia=1:numel(M_final)
    for ii=1:numel(M_final{ia})
        ColorAllNeuronsForMo(M_final{ia}{ii},300,400,['Day' num2str(ii) 'after Day' num2str(ia)], 'X:\EmilyShijieShared_old\6922_moBatchVer',xxsfyysf,[10,10]);
    %ColorAllNeuronsForMo(A2,300,400,['18' 'Newupdating after' num2str(ia)],outputdirForA,xxsfyysf);
    end
end

M=M_final;
%% Save results
Vars = 'M_final shifts template birdnum';
eval(sprintf('save %scnmfe_BatchVer_PartII_MotionCorrection.mat %s -v7.3', cnmfedir, Vars));
