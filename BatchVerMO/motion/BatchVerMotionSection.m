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

%[~, ~, ~] = maybe_spawn_workers(4); 

%%
global shifts
M=cell(1,AnumInFolder);  M1=cell(1,AnumInFolder); 
shifts_consec=cell(1,AnumInFolder); shifts_consec_up_1=cell(1,AnumInFolder); shifts_consec_up_2=cell(1,AnumInFolder);
shifts_inter=cell(1,AnumInFolder);  shifts_inter_up_1=cell(1,AnumInFolder);  shifts_inter_up_2=cell(1,AnumInFolder);
ind_del=cell(1,AnumInFolder);

for ia=2:AnumInFolder
   
    Y_previousday=Y(:,:,ia-1);%Y_ex_oneday=Y(:,:,reseq_ind(2:end));
    Y_oneday=Y(:,:,ia);
    As_oneday=AsfromDaysCell{ia}; %As_ex_oneday=AsfromDaysCell(reseq_ind(2:end));
    siz_oneday=sizes(ia);         %siz_ex_oneday=sizes(reseq_ind(2:end));
    
    gridstartend=[1,300,61,400,1,1];
    %tic; [M{ia},shifts{ia},~,xxsfyysf,ind_del{ia}] = normcorre_BatchVer(Y_ex_oneday,options_nonrigid,Y_oneday,siz_ex_oneday,As_ex_oneday,startendgrid,update_num); toc
    tic; [M_consec_temp,shifts_consec{ia},shifts_consec_up_1{ia},shifts_consec_up_2{ia},~,xxsfyysf,ind_del_temp] = ...
        normcorre_BatchVer(Y_oneday,options_nonrigid,Y_previousday,siz_oneday,As_oneday,gridstartend,update_num); toc
    M{ia-1}{ia}=M_consec_temp{1};
    ind_del{ia-1}{ia}=ind_del_temp{1};

    % Based on consecutive days' shift matrix, fill in all motions
    % The day's own data
    M{ia}{ia}=As_oneday;
    ind_del{ia}{ia}=false(1,siz_oneday);
    
    % inverse consec
    Mf_temp=[];
    ind_del{ia}{ia-1}=false(1,sizes(ia-1));
    for ni=1:sizes(ia-1)
        Y_one_neuron=reshape(AsfromDaysCell{ia-1}(:,ni),size(Y,1),size(Y,2));
        Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)) = ...
            imwarp(Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)),-cat(3,shifts_consec_up_2{ia}.*(-1),shifts_consec_up_1{ia}.*(-1)),options_nonrigid.shifts_method);
        if any(Y_one_neuron)==0
            ind_del{ia}{ia-1}(ni)=true;
        end
        Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
    end
    M{ia}{ia-1}=Mf_temp;
    
         
    if ia>=3
        % Update data, fill the registration to the previous 2nd day
        Y_previous2day=Y(:,:,ia-2);
        currentday=A2image(M{ia-1}{ia});
        Y_2day=mean(Y(:,:,ia-1),currentday,3);
        
        
        tic; [M_inter_temp,shifts_inter{ia},shifts_inter_up_1{ia},shifts_inter_up_2{ia},~,xxsfyysf,ind_del_temp] = ...
            normcorre_BatchVer(Y_2day,options_nonrigid,Y_previous2day,siz_oneday,M{ia-1}{ia},gridstartend,update_num); toc
        M{ia-2}{ia}=M_inter_temp{1};
        ind_del{ia-2}{ia}=ind_del_temp{1};
        
        % inverse inter
        Mf_temp=[];
        ind_del{ia}{ia-2}=false(1,sizes(ia-2));
        for ni=1:sizes(ia-2)
            Y_one_neuron=reshape(AsfromDaysCell{ia-2}(:,ni),size(Y,1),size(Y,2));
            Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)) = ...
                imwarp(Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)),-cat(3,shifts_inter_up_2{ia}.*(-1),shifts_inter_up_1{ia}.*(-1)),options_nonrigid.shifts_method);
            if any(Y_one_neuron)==0
                ind_del{ia}{ia-2}(ni)=true;
            end
            Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
        end
        M{ia}{ia-2}=Mf_temp;
    end
    
    if ia>=4
        shifts_1=cumsum(flip(cat(3,shifts_up_1{2:ia}),3),3);
        shifts_2=cumsum(flip(cat(3,shifts_up_2{2:ia}),3),3);
        % for previous days, adding a new day means adding a new day to register 
        % and for current new day, previous days need to be registered.
        for io=1:ia-3
            % for previous days
            Mf_temp=[];
            ind_del{io}{ia} = false(1,siz_oneday);
            shifts_1_temp=shifts_1(:,:,ia-io);
            shifts_2_temp=shifts_2(:,:,ia-io);
            for ni=1:siz_oneday
                Y_one_neuron=reshape(As_oneday(:,ni),size(Y,1),size(Y,2));
                Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)) = imwarp(Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)),-cat(3,shifts_2_temp,shifts_1_temp),options_nonrigid.shifts_method);
                if any(Y_one_neuron)==0
                    ind_del{io}{ia}(ni)=true;
                end
                Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
            end
            M{io}{ia}=Mf_temp;
            
            % for current day
            Mf_temp=[];
            ind_del{ia}{io} = false(1,sizes(io));
            for ni=1:sizes(io)
                Y_one_neuron=reshape(AsfromDaysCell{io}(:,ni),size(Y,1),size(Y,2));
                Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)) = imwarp(Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)),-cat(3,shifts_2_temp.*(-1),shifts_1_temp.*(-1)),options_nonrigid.shifts_method);
                if any(Y_one_neuron)==0
                    ind_del{ia}{io}(ni)=true;
                end
                Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
            end
            M{ia}{io}=Mf_temp;
        end
    end
end
% the only left out one;
M{1}{1}=AsfromDaysCell{1};  ind_del{1}{1}=false(1,sizes(1));  
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
