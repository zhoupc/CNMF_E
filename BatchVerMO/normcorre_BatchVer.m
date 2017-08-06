function [M_final,shifts_g,shifts_g_up_1,shifts_g_up_2,template,xxsfyysf,ind_del_full] = normcorre_BatchVer(Y,options,template,sizes,As,gridstartend,update_num)

% online motion correction through DFT subpixel registration
% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction (optional, rigid registration is performed if not provided)
% template:         provide template (optional)
% sizes:            cell array, sizes{i}=size(Ai,2);

% OUTPUTS
% M_final:          motion corrected data
% shifts_up:        upsampled shifts
% shifts:           originally calculated shifts
% template:         calculated template

%% first determine filetype

if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff');
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        sizY = size(Y.Y);
        T = sizY(end);
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5');
        filetype = 'hdf5';
        fileinfo = hdf5info(Y);
        data_name = fileinfo.GroupHierarchy.Datasets.Name;
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        T = sizY(end);
    elseif strcmpi(ext,'raw')
        filetype = 'raw';
        fid = fopen(Y);
        FOV = [options.d1,options.d2];
        bitsize = options.bitsize;
        imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        T = file_length/imsize;
        sizY = [FOV,T];
        fclose(fid);        
    end    
elseif isobject(Y);
    filetype = 'mem';
    sizY = size(Y,'Y');
    T = sizY(end);
else % array loaded in memory
    filetype = 'mat';
    %Y = double(Y);
    sizY = size(Y);
    T = sizY(end);
end

nd = length(sizY)-1;                          % determine whether imaging is 2d or 3d
if or(nd==2,nd==3)
    sizY = sizY(1:nd);
elseif nd==1
    T=1;
end
%% set default parameters if not present

if ~exist('options','var') || isempty(options);
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2));
    if nd > 2; options.d3 = sizY(3); end
end

memmap = options.memmap;
grid_size = options.grid_size; 
mot_uf = options.mot_uf;
min_patch_size = options.min_patch_size;
overlap_pre = options.overlap_pre;
overlap_post = options.overlap_post;
upd_template = options.upd_template;
bin_width = options.bin_width;
buffer_width = options.buffer_width;
max_dev_g = options.max_dev;
init_batch = options.init_batch;
us_fac = options.us_fac;
method = options.method;
filename = options.mem_filename;
iter = options.iter;
add_value = options.add_value;
max_shift = options.max_shift;

while mod(T,bin_width) == 1
    if T == 1
        error('Movie appears to have only one frame. Use the function normcorre instead')        
    end
    bin_width = bin_width + 1;
end

%% first check for offset due to bi-directional scanning

if options.correct_bidir
    col_shift = correct_bidirectional_offset(Y,options.nFrames,options.bidir_us);
    if col_shift
        if strcmpi(options.shifts_method,'fft')
            options.shifts_method = 'cubic';
            fprintf('Offset %1.1f pixels due to bidirectional scanning detected. Cubic shifts will be applied. \n',col_shift); 
        end
    end
else
    col_shift = 0;
end


%% read initial batch and compute template

init_batch = min(T,init_batch);
perm = randperm(T,init_batch);
switch filetype
    case 'tif'
        Y1 = imread(Y,'Index',perm(1),'Info',tiffInfo);
        Y_temp = zeros(sizY(1),sizY(2),init_batch,'like',Y1);
        Y_temp(:,:,1) = Y1;
        for tt = 2:init_batch
            Y_temp(:,:,tt) = imread(Y,'Index',perm(tt),'Info',tiffInfo);
        end
    case 'hdf5'
        Y_temp = bigread2(Y,1,init_batch);        
    case 'mem'
        if nd == 2; Y_temp = Y.Y(:,:,1:init_batch); elseif nd == 3; Y_temp = Y.Y(:,:,:,1:init_batch); end
    case 'mat'
        if nd == 2; Y_temp = Y(:,:,perm); elseif nd == 3; Y_temp = Y(:,:,:,perm);
        elseif nd==1; Y_temp = Y(:,:); end
    case 'raw'
         Y_temp = read_raw_file(Y,1,init_batch,FOV,bitsize);
end
data_type = class(Y_temp);
Y_temp = single(Y_temp);

if nargin < 3 || isempty(template)
    template_in = median(Y_temp,nd+1)+add_value;
else
    template_in = single(template + add_value);
end

[d1,d2,d3,~] = size(Y_temp);
if or(nd==1,nd == 2); d3 = 1; end
%% setup grids for patches

[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size,gridstartend);
shifts_g = struct('shifts',cell(T,1),'shifts_up',cell(T,1),'diff',cell(T,1));
temp_cell = mat2cell_ov(template_in,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,2);
xxsfyysf{1}=xx_s;   xxsfyysf{2}=xx_f;   xxsfyysf{3}=yy_s;   xxsfyysf{4}=yy_f;
%% precompute some quantities that are used repetitively for template matching and applying shifts
Nr = cell(size(temp_cell));
Nc = cell(size(temp_cell));
Np = cell(size(temp_cell));
Bs = cell(size(temp_cell));
for i = 1:length(xx_us)
    for j = 1:length(yy_us)
        for k = 1:length(zz_us)
            [nr,nc,np] = size(temp_cell{i,j,k});
            nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
            nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
            np = ifftshift(-fix(np/2):ceil(np/2)-1);
            [Nc{i,j,k},Nr{i,j,k},Np{i,j,k}] = meshgrid(nc,nr,np);
            extended_grid = [max(xx_us(i)-overlap_post(1),gridstartend(1)),min(xx_uf(i)+overlap_post(1),gridstartend(2)),max(yy_us(j)-overlap_post(2),gridstartend(3)),min(yy_uf(j)+overlap_post(2),gridstartend(4)),max(zz_us(k)-overlap_post(3),gridstartend(5)),min(zz_uf(k)+overlap_post(3),gridstartend(6))];            
            Bs{i,j,k} = permute(construct_weights([xx_us(i),xx_uf(i),yy_us(j),yy_uf(j),zz_us(k),zz_uf(k)],extended_grid),[2,1,3]); 
        end
    end
end
if or(nd == 2,nd==1); Np = cellfun(@(x) 0,Nr,'un',0); end

%%
maxNumCompThreads(1);
template = mat2cell_ov(template_in,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,2);
temp_mat = template_in(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6));
display('Originalline176')
display(size(temp_mat))
use_windowing = options.use_windowing;
phase_flag = options.phase_flag;

if use_windowing
    fftTemp = cellfun(@fftn,cellfun(@han,template,'un',0),'un',0);
    fftTempMat = fftn(han(temp_mat));
else
    fftTemp = cellfun(@fftn,template,'un',0);
    fftTempMat = fftn(temp_mat);
end

if nd==1||nd == 2; buffer = mat2cell_ov(zeros(d1,d2,bin_width),xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,2); end
if nd == 3; buffer = mat2cell_ov(zeros(d1,d2,d3,bin_width),xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,3); end

if ~strcmpi(options.output_type,'mat')
    options.mem_batch_size = max(min(round(options.mem_batch_size/bin_width)*bin_width,T),1);
    if nd == 2; mem_buffer = zeros(d1,d2,options.mem_batch_size,'single'); end
    if nd == 3; mem_buffer = zeros(d1,d2,d3,options.mem_batch_size,'single'); end
end

switch lower(options.output_type)
    case 'mat'
        %M_final = zeros([sizY,T],data_type);
        M_final = cell(1,T);
    case 'memmap'
        M_final = matfile(filename,'Writable',true);
        if nd == 2; M_final.Y(d1,d2,T) = zeros(1,data_type); end
        if nd == 3; M_final.Y(d1,d2,d3,T) = zeros(1,data_type); end
        M_final.Yr(d1*d2*d3,T) = zeros(1,data_type);        
    case {'hdf5','h5'}
         if exist(options.h5_filename,'file')
            [pathstr,fname,ext] = fileparts(options.h5_filename);             
            new_filename = fullfile(pathstr,[fname,'_',datestr(now,30),ext]);
            warning_msg = ['File ',options.h5_filename,'already exists. Saving motion corrected file as',new_filename];            
            warning('%s',warning_msg);
            options.h5_filename = new_filename;
        end       
        M_final = options.h5_filename;
        if nd == 2
            h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,options.mem_batch_size],'Datatype',data_type);
        elseif nd == 3
            h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,d3,Inf],'Chunksize',[d1,d2,d3,options.mem_batch_size],'Datatype',data_type);
        end
    case {'tif','tiff'}
        M_final = ['motion corrected file has been saved as ', options.tiff_filename];
        opts_tiff.append = true;
        opts_tiff.big = true;
        if nd == 3
            error('Saving volumetric tiff stacks is currently not supported. Use a different filetype');
        end
    otherwise
        error('This filetype is currently not supported')
end   

cnt_buf = 0;
fprintf('Template initialization complete. \n')
%%
ind_del_full = cell(1,T);
for it = 1:iter
    for t = 1:bin_width:T
        tpointer=t;
        switch filetype
            case 'tif'
                Ytm = zeros(sizY(1),sizY(2),min(t+bin_width-1,T)-t+1,'single');
                for tt = 1:min(t+bin_width-1,T)-t+1
                    Ytm(:,:,tt) = single(imread(Y,'Index',t+tt-1,'Info',tiffInfo));
                end
            case 'hdf5'
                Ytm = single(h5read(Y,data_name,[ones(1,nd),t],[sizY(1:nd),min(t+bin_width-1,T)-t+1]));
            case 'mem'
                if nd == 2; Ytm = single(Y.Y(:,:,t:min(t+bin_width-1,T))); end
                if nd == 3; Ytm = single(Y.Y(:,:,:,t:min(t+bin_width-1,T))); end
            case 'mat'
                if nd == 1; Ytm = single(Y(:,:)); end
                if nd == 2; Ytm = single(Y(:,:,t:min(t+bin_width-1,T))); end
                if nd == 3; Ytm = single(Y(:,:,:,t:min(t+bin_width-1,T))); end
            case 'raw'
                Ytm = single(read_raw_file(Y,t,min(t+bin_width-1,T)-t+1,FOV,bitsize));                
        end
        
        if nd == 2||nd == 1; 
            nd=2;
            if size(Ytm,3)==1
                Ytc = mat2cell(Ytm,d1,d2);
            else
                Ytc = mat2cell(Ytm,d1,d2,ones(1,size(Ytm,ndims(Ytm))));    end
        elseif nd == 3; 
            if size(Ytm,4)==1
                Ytc = mat2cell(Ytm,d1,d2,d3);
            else
                Ytc = mat2cell(Ytm,d1,d2,d3,ones(1,size(Ytm,ndims(Ytm)))); end
        end
        
        Mf = cell(size(Ytc));
        lY = length(Ytc);
        display('lY=')
        display(lY)
        ind_del = cell(1,lY);
        %buffer = cell(length(xx_us),length(yy_us),length(zz_us),size(Ytm,ndims(Ytm)));
        shifts = struct('shifts',cell(lY,1),'shifts_up',cell(lY,1),'diff',cell(lY,1));
        %buf = struct('Mf',cell(lY,1));
        parfor ii = 1:lY
            display('line 281')
            Yt = Ytc{ii}; 
            Yt_trunked = Yt(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6));
            minY = min(Yt(:));
            maxY = max(Yt(:));
            Yc = mat2cell_ov(Yt,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,2);
            
            if use_windowing
                fftY = cellfun(@fftn, cellfun(@han,Yc, 'un',0),'un',0);
            else
                fftY = cellfun(@fftn, Yc, 'un',0);
            end
            
            neuronnumber = sizes(ii+tpointer-1);
                                         %M_fin = cell(length(sizes(ii)),length(xx_us),length(yy_us),length(zz_us)); %zeros(size(Y_temp));
            M_fin = cell(1,neuronnumber);%struct('M_fin',cell(length(xx_us),length(yy_us),length(zz_us)));
            shifts_temp = zeros(length(xx_s),length(yy_s),length(zz_s),nd); 
            diff_temp = zeros(length(xx_s),length(yy_s),length(zz_s));
            display('line299')
            if numel(shifts_temp) > 1      
                if use_windowing
                    if nd == 2; out_rig = dftregistration_min_max(fftTempMat,fftn(han(Yt_trunked)),us_fac,-max_shift,max_shift,phase_flag); lb = out_rig(3:4); ub = out_rig(3:4); end
                    if nd == 3; out_rig = dftregistration_min_max_3d(fftTempMat,fftn(han(Yt_trunked)),us_fac,-max_shift,max_shift,phase_flag); lb = out_rig(3:5); ub = out_rig(3:5); end
                else
                    if nd == 2; out_rig = dftregistration_min_max(fftTempMat,fftn(Yt_trunked),us_fac,-max_shift,max_shift,phase_flag); lb = out_rig(3:4); ub = out_rig(3:4); end
                    if nd == 3; out_rig = dftregistration_min_max_3d(fftTempMat,fftn(Yt_trunked),us_fac,-max_shift,max_shift,phase_flag); lb = out_rig(3:5); ub = out_rig(3:5); end
                end
                max_dev = max_dev_g;
            else
                lb = -max_shift(1,nd);
                ub = max_shift(1,nd);
                max_dev = 0*max_dev_g;
            end
            for i = 1:length(xx_s)
                for j = 1:length(yy_s)           
                    for k = 1:length(zz_s)
                        if nd == 2
                            %[output,Greg] = dftregistration_min_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev(1:2),ub+max_dev(1:2));  
                            output = dftregistration_min_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev(1:2),ub+max_dev(1:2),phase_flag);  
                        elseif nd == 3
                            output = dftregistration_min_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev,ub+max_dev,phase_flag); 
                            shifts_temp(i,j,k,3) = output(5);
                        end
                       
                        shifts_temp(i,j,k,1) = output(3);
                        shifts_temp(i,j,k,2) = output(4); 
                        diff_temp(i,j,k) = output(2);
                        if all([length(xx_s),length(yy_s),length(zz_s)] == 1) && strcmpi(options.shifts_method,'fft')
                            for ni=1:sizes(ii+tpointer-1)
                                Y_one_neuron=reshape(As(:,ni),d1,d2);
                                Y_one_neuron_c = mat2cell_ov(Y_one_neuron,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,2);
                                M_fin{ni}{i,j,k} = shift_reconstruct(Y_one_neuron_c{i,j,k},shifts_temp(i,j,k,:),diff_temp(i,j,k),us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k},options.boundary,add_value);
                            end
                        end                                               
                    end
                end
            end            
            
            shifts(ii).shifts = shifts_temp;
            shifts(ii).diff = diff_temp;
            display('line340')
            switch lower(options.shifts_method)
                case 'fft'
                    display('Using fft.')
                    if any([length(xx_s),length(yy_s),length(zz_s)] > 1)
                        if mot_uf(3) > 1
                            tform = affine3d(diag([mot_uf(:);1]));
                            diff_up = imwarp(diff_temp,tform,'OutputView',imref3d([length(xx_uf),length(yy_uf),length(zz_uf)]));
                            shifts_up = zeros([size(diff_up),3]);
                            for dm = 1:3; shifts_up(:,:,:,dm) = imwarp(shifts_temp(:,:,:,dm),tform,'OutputView',imref3d([length(xx_uf),length(yy_uf),length(zz_uf)])); end
                        else
                            shifts_up = imresize(shifts_temp,[length(xx_uf),length(yy_uf)]);
                            diff_up = imresize(diff_temp,[length(xx_uf),length(yy_uf)]);
                        end
                        shifts(ii).shifts_up = shifts_up;
                        shifts(ii).diff = diff_up;
                        for ni=1:sizes(ii+tpointer-1)
                            for i = 1:length(xx_uf)
                                for j = 1:length(yy_uf)
                                    for k = 1:length(zz_uf)
                                        Y_one_neuron=reshape(As(:,ni),d1,d2);
                                        extended_grid_2 = [max(xx_us(i)-overlap_post(1),gridstartend(1)),min(xx_uf(i)+overlap_post(1),gridstartend(2)),max(yy_us(j)-overlap_post(2),gridstartend(3)),min(yy_uf(j)+overlap_post(2),gridstartend(4)),max(zz_us(k)-overlap_post(3),gridstartend(5)),min(zz_uf(k)+overlap_post(3),gridstartend(6))];
                                        I_temp = Y_one_neuron(extended_grid_2(1):extended_grid_2(2),extended_grid_2(3):extended_grid_2(4),extended_grid_2(5):extended_grid_2(6));
                                        %I_temp = Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6));
                                        M_fin{ni}{i,j,k} = shift_reconstruct(I_temp,shifts_up(i,j,k,:),diff_up(i,j,k),us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k},options.boundary,add_value);
                                        %M_fin{i,j,k} = shift_reconstruct2(I_temp,shifts_up(i,j,k,:),'bilinear',diff_up(i,j,k),us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k},options.boundary,add_value);
                                    end
                                end
                            end
                        end
                    else
                        shifts_up = shifts_temp;
                        shifts(ii).shifts_up = shifts(ii).shifts;
                    end
                    gx = max(abs(reshape(diff(shifts_up,[],1),[],1)));
                    gy = max(abs(reshape(diff(shifts_up,[],2),[],1)));
                    gz = max(abs(reshape(diff(shifts_up,[],3),[],1)));
                    flag_interp = max([gx;gy;gz;0])<0.5;      % detect possible smearing
                    
                    %Mf_each=cell(1,sizes(ii));
                    Mf_temp=[];
                    ind_del{ii} = false(1,sizes(ii+tpointer-1));
                    for ni=1:sizes(ii+tpointer-1)
                        M_fin_tmp=M_fin{ni};
                        if flag_interp
                            Mf_each = cell2mat_ov_sum(M_fin_tmp,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,gridstartend,Bs) - add_value;
                        else
                            Mf_each = cell2mat_ov(M_fin_tmp,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,gridstartend) - add_value;
                        end
                        Mf_each(Mf_each<minY)=minY;
                        Mf_each(Mf_each>maxY)=maxY;
                        Y_one_neuron=reshape(As(:,ni),d1,d2);
                        Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6))=Mf_each;
                        if any(Y_one_neuron)==0
                            ind_del{ii}(ni)=true;
                        end
                        Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
                    end
                    Mf{ii}=Mf_temp;
                otherwise
                    display('Shift using other methods.')
                    %shifts(ii).shifts_up = shifts(ii).shifts;
                    if nd == 3
                        tform = affine3d(diag([mot_uf(:);1]));
                        shifts_up = zeros([options.d1,options.d2,options.d3,3]);
                        for dm = 1:3; shifts_up(:,:,:,dm) = imwarp(shifts_temp(:,:,:,dm),tform,'OutputView',imref3d([options.d1,options.d2,options.d3])); end
                        shifts_up(2:2:end,:,:,2) = shifts_up(2:2:end,:,:,2) + col_shift;
                        Mf{ii} = imwarp(Yt,-cat(3,shifts_up(:,:,2),shifts_up(:,:,1)),options.shifts_method);
                    else
                        display('line410')
                        shifts_up = imresize(shifts_temp,[gridstartend(2)-gridstartend(1)+1,gridstartend(4)-gridstartend(3)+1]);
                        %shifts(ii).shifts_up = shifts_up;
                        shifts_up(2:2:end,:,2) = shifts_up(2:2:end,:,2) + col_shift;
                        shifts(ii).shifts_up = shifts_up;
                        Mf_temp=[];
                        display('line416')
                        ind_del{ii} = false(1,sizes(ii+tpointer-1));
                        for ni=1:sizes(ii+tpointer-1)
                            Y_one_neuron=reshape(As(:,ni),d1,d2);
                            Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)) = imwarp(Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)),-cat(3,shifts_up(:,:,2),shifts_up(:,:,1)),options.shifts_method);
                            if any(Y_one_neuron)==0
                                ind_del{ii}(ni)=true;
                            end
                            Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
                        end
                        display('line424')
                        Mf{ii}=Mf_temp;
                    end
            end
        end

        shifts_g(t:min(t+bin_width-1,T)) = shifts;
        %shifts_g_up_1=zeros(size(shifts(1).shifts_up,1),size(shifts(1).shifts_up,2),T); shifts_g_up_2=zeros(size(shifts(1).shifts_up,1),size(shifts(1).shifts_up,2),T);
        shifts_g_up_1(:,:) = shifts(1).shifts_up(:,:,1);
        shifts_g_up_2(:,:) = shifts(1).shifts_up(:,:,2);
        %Mf = cell2mat(Mf);
        %%%%%%%%%%% Updating datasets%%%%%%%%%%%%
        MfMAT=[];
        for ii=1:numel(Mf)
            k = size(Mf{ii},2);  C=ones(k,1);
            Brainbow = Mf{ii}*C; Brainbow = reshape(Brainbow,300,400);
            MfMAT=cat(3,MfMAT,Brainbow);
        end
        %Y(:,:,t:min(t+bin_width-1,T))=MfMAT;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~strcmpi(options.output_type,'mat')
            rem_mem = rem(t+lY-1,options.mem_batch_size);
            if rem_mem == 0; rem_mem = options.mem_batch_size; end            
            if nd == 2; mem_buffer(:,:,rem_mem-lY+1:rem_mem) = cast(Mf,data_type); end
            if nd == 3; mem_buffer(:,:,:,rem_mem-lY+1:rem_mem) = cast(Mf,data_type); end
        end
        if it == iter
        switch lower(options.output_type)
            case 'mat'
                if nd == 2; M_final(t:min(t+bin_width-1,T)) = Mf;  
                            ind_del_full(t:min(t+bin_width-1,T)) = ind_del; end %cast(Mf,data_type); end
                if nd == 3; M_final(:,:,:,t:min(t+bin_width-1,T)) = cast(Mf,data_type); end
            case 'memmap'
                if rem_mem == options.mem_batch_size || t+lY-1 == T
                    if nd == 2; M_final.Y(:,:,t+lY-rem_mem:t+lY-1) = mem_buffer(:,:,1:rem_mem); end
                    if nd == 3; M_final.Y(:,:,:,t+lY-rem_mem:t+lY-1) = mem_buffer(:,:,:,1:rem_mem); end
                    M_final.Yr(:,t+lY-rem_mem:t+lY-1) = reshape(mem_buffer(1:d1*d2*d3*rem_mem),d1*d2*d3,rem_mem);
                end      
            case {'hdf5','h5'}
                if rem_mem == options.mem_batch_size || t+lY-1 == T
                    if nd == 2; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,1:rem_mem),[ones(1,nd),t+lY-rem_mem],[sizY(1:nd),rem_mem]); end
                    if nd == 3; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,:,1:rem_mem),[ones(1,nd),t+lY-rem_mem],[sizY(1:nd),rem_mem]); end
                end
            case {'tif','tiff'}
                if rem_mem == options.mem_batch_size || t+lY-1 == T
                    saveastiff(cast(mem_buffer(:,:,1:rem_mem),data_type),options.tiff_filename,opts_tiff);
                end
        end        
        end
        
        % update template
        fprintf('%i out of %i frames registered, iteration %i out of %i \n',t+lY-1,T,it,iter)
        if and(upd_template,t<=update_num)
            display('Here')
            display(t)
            display('Updating Template Here')
            cnt_buf = cnt_buf + 1;
            [~,~,~,~,zz_s_temp,zz_f_temp,~,~,~,~,~,~] = construct_grid([grid_size(1) grid_size(2) size(MfMAT,3)],mot_uf,d1,d2,size(MfMAT,3), min_patch_size,[gridstartend(1:5), size(MfMAT,3)]);
            buffer = mat2cell_ov(MfMAT,xx_s,xx_f,yy_s,yy_f,zz_s_temp,zz_f_temp,overlap_pre,3);
            if strcmpi(method{2},'mean')
                new_temp = cellfun(@(x) nanmean(x,nd+1), buffer, 'UniformOutput',false);
            elseif strcmpi(method{2},'median');
                new_temp = cellfun(@(x) nanmedian(x,nd+1), buffer, 'UniformOutput', false);
            end
            if any(reshape(cell2mat(cellfun(@(x) any(isnan(x(:))), new_temp, 'un',0)),[],1))
                parfor i = 1:numel(new_temp)
                    new_temp{i}(isnan(new_temp{i})) =  template{i}(isnan(new_temp{i}));
                end
            end
            if strcmpi(method{1},'mean')
                cnt = t/bin_width + 1;
                template = cellfun(@plus, cellfun(@(x) x*(cnt-1)/cnt, template,'un',0), cellfun(@(x) x*1/cnt, new_temp,'un',0), 'un',0);
            elseif strcmpi(method{1},'median');
                if cnt_buf <= buffer_width
                    if nd == 2; buffer_med(:,:,cnt_buf) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,gridstartend); end
                    if nd == 3; buffer_med(:,:,:,cnt_buf) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,gridstartend); end
                else
                    buffer_med = circshift(buffer_med,[zeros(1,nd),-1]);
                    if nd == 2; buffer_med(:,:,buffer_width) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,gridstartend); end
                    if nd == 3; buffer_med(:,:,:,buffer_width) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,gridstartend); end
                end
                template = mat2cell_ov(nanmedian(buffer_med,nd+1),xx_s-gridstartend(1)+1,xx_f-gridstartend(1)+1,yy_s-gridstartend(3)+1,yy_f-gridstartend(3)+1,zz_s-gridstartend(5)+1,zz_f-gridstartend(5)+1,overlap_pre,2);
            end
            temp_mat = cell2mat_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,gridstartend);
            
            if use_windowing
                fftTemp = cellfun(@fftn, cellfun(@han,template, 'un',0),'un',0);            
                fftTempMat = fftn(han(temp_mat));            
            else
                fftTemp = cellfun(@fftn, template, 'un',0);            
                fftTempMat = fftn(temp_mat);
            end
            display(size(fftTempMat))
        end
        
    end

    if it == iter
        template = cellfun(@(x) x - add_value,template,'un',0);
        template = cell2mat_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,gridstartend);
        
    end
    if memmap;
        M_final.shifts = shifts_g;
        M_final.template = template;
    end
    maxNumCompThreads(1);
end