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
AsfromDaysCell=cell(1,AnumInFolder);
AsfromDaysCell_central=cell(1,AnumInFolder);
sizes=[];
for ia=1:AnumInFolder
    Anext=load(fullfile(cnmfedir,Alist(ia).name));
    %ColorAllNeurons(Anext.Afinal,300,400,num2str(ia),outputdirForA)
    Atemp=Anext.Afinal;
    AsfromDaysCell{ia}=Atemp;   % individual column is each A, for actual motion correction

    Atemp=centralA(Atemp);
    AsfromDaysCell_central{ia}=Atemp;

    k = size(Atemp,2);
    sizes=[sizes k];
    AsfromDaysPic=cat(3,AsfromDaysPic,A2image(Atemp,300,400)); % whole picture, for registering, getting shifts        
end

%% 2. for each A in AsfromDays, register others to it.
Y = AsfromDaysPic;
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));

options_nonrigid = NoRMCorreSetParms('upd_template',false,'iter',1,...
                                     'd1',size(Y,1),'d2',size(Y,2),'grid_size',[80,80],'min_patch_size',[50,50,1],'overlap_pre',[10,10,1],'overlap_post',[10,10,1],...
                                     'mot_uf',10,'bin_width',1,...
                                     'shifts_method','cubic',...
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
gridstartend=[1,300,61,400,1,1];
%%
M=cell(1,AnumInFolder);  M_central=cell(1,AnumInFolder); 
ind_del=cell(1,AnumInFolder);
% the only left out one;
M{1}{1}=AsfromDaysCell{1};  ind_del{1}{1}=false(1,sizes(1));
M_central{1}{1}=centralA(AsfromDaysCell{1});

for ia=2:AnumInFolder
    siz_oneday=sizes(ia);
    % The day's own data
    M{ia}{ia}=AsfromDaysCell{ia};
    M_central{ia}{ia}=centralA(AsfromDaysCell{ia});
    ind_del{ia}{ia}=false(1,siz_oneday);
    % Previous days
    M_buffer=[];
    for io=ia-1:-1:1
        if io==ia-1
            Y_template=Y(:,:,io);
        else
            Y_template=A2image(cat(2,M{io}{io:ia-1}),size(Y,1),size(Y,2));
        end
        if io==ia-1
            Y_toregister=Y(:,:,ia);
        else
            Y_toregister=A2image(cat(2,M{io+1}{io+1:ia}),size(Y,1),size(Y,2));
        end
 
        %tic; [M{ia},shifts{ia},~,xxsfyysf,ind_del{ia}] = normcorre_BatchVer(Y_ex_oneday,options_nonrigid,Y_oneday,siz_ex_oneday,As_ex_oneday,startendgrid,update_num); toc
        tic; [M_temp,shifts,shifts_up_1,shifts_up_2,~,xxsfyysf,ind_del_temp] = ...
            normcorre_BatchVer(Y_toregister,options_nonrigid,Y_template,siz_oneday,M{io+1}{ia},gridstartend,update_num); toc
        M{io}{ia}=M_temp{1};
        M_central{io}{ia}=centralA(M_temp{1});
        ind_del{io}{ia}=ind_del_temp{1};

        % inverse consec
        Mf_temp=[];
        ind_del{ia}{io}=false(1,sizes(io));
        for ni=1:sizes(io)
            Y_one_neuron=reshape(AsfromDaysCell{io}(:,ni),size(Y,1),size(Y,2));
            Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)) = ...
                imwarp(Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)),-cat(3,shifts_up_2.*(-1),shifts_up_1.*(-1)),options_nonrigid.shifts_method);
            if any(Y_one_neuron)==0
                ind_del{ia}{io}(ni)=true;
            end
            Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
        end
        M{ia}{io}=Mf_temp;
        M_central{ia}{io}=centralA(M{ia}{io});
    end
end

%%
ind_del_full_cell=cellfun(@(x) cell2mat(x), ind_del, 'UniformOutput',0);
ind_del_full=sum(reshape(cell2mat(ind_del_full_cell),1,sum(sizes),[]),3)>0;    %sum(cat(3,ind_del{:}))>0;
ind_del_full_cell=mat2cell(~ind_del_full,1,sizes);
N_eachday=cellfun(@(x) sum(x), ind_del_full_cell);
M_concat=cellfun(@(x) cat(2,x{:}), M, 'UniformOutput',0);
M_del = cellfun(@(x) x(:,~ind_del_full), M_concat, 'UniformOutput',0);
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
