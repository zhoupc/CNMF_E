function [imData,sizx,sizy]=smod_bigread2(varargin)
%reads tiff files in Matlab bigger than 4GB, allows reading from sframe to sframe+num2read-1 frames of the tiff - in other words, you can read page 200-300 without rading in from page 1.
%based on a partial solution posted on Matlab Central (http://www.mathworks.com/matlabcentral/answers/108021-matlab-only-opens-first-frame-of-multi-page-tiff-stack)
%Darcy Peterka 2014, v1.0
%Darcy Peterka 2014, v1.1
%Darcy Peterka 2016, v1.2
%Eftychios Pnevmatikakis 2016, v1.3 (added hdf5 support)
%Darcy Peterka 2016, v1.4
%Darcy Peterka 2016, v1.5(bugs to dp2403@columbia.edu)
%Program checks for bit depth, whether int or float, and byte order.  Assumes uncompressed, non-negative (i.e. unsigned) data.
%
% Usage:  my_data=bigread('path_to_data_file, start frame, num to read,big);
% "my_data" will be your [M,N,frames] array.
%Will do mild error checking on the inputs - last three inputs are optional -
%if nargin == 2, then assumes second is number of frames to read, and reads
%till the end
%Checks to see if imageJ generated tif, then uses info in imageJ generated
%image description to load the files based on offset and image number.
%

%modified by Pengcheng Zhou, Colubmia University, 2017
% it supports following data formats:
% tiff
% hdf5
% avi

%get image info
path_to_file=strjoin(varargin(1));

if nargin<2
    sframe = 1;
else
    sframe = max(1, round(varargin{2}));
end
[~,~,ext] = fileparts(path_to_file);

if strcmpi(ext,'.tiff') || strcmpi(ext,'.tif') || strcmpi(ext,'.BTF');
    info = imfinfo(path_to_file);
    
    
    if isfield(info(1), 'ImageDescription') && ~isempty(strfind(info(1).ImageDescription,'ImageJ'))
        junk1=regexp(info(1).ImageDescription,'images=\d*','match');
        if isempty(junk1)
            junk1=regexp(info(1).ImageDescription,'frames=\d*','match');
        end
        junk2=strjoin(junk1);
        aa=textscan(junk2,'%*s %d','delimiter','=');
        aa = aa{1};
        
        numFrames=aa;
        num_tot_frames=numFrames;
        
        if nargin<3
            if nargin<2
                sframe = 1;
            else
                sframe=cell2mat(varargin(2));
            end
            num2read=numFrames-sframe+1;
        end
        if nargin==3
            num2read=cell2mat(varargin(3));
            sframe = cell2mat(varargin(2));
        end
        if sframe<=0
            sframe=1;
        end
        if num2read<1
            num2read=1;
        end
        %sframe = 1;
        %num2read=numFrames-sframe+1;
        
        if sframe>num_tot_frames
            sframe=num_tot_frames;
            num2read=1;
            display('starting frame has to be less than number of total frames...');
        end
        if (num2read+sframe<= num_tot_frames+1)
            lastframe=num2read+sframe-1;
        else
            num2read=numFrames-sframe+1;
            lastframe=num_tot_frames;
            display('Hmmm...just reading from starting frame until the end');
        end
        
        
        bd=info.BitDepth;
        he=info.ByteOrder;
        bo=strcmp(he,'big-endian');
        if (bd==64)
            form='double';
            for_mult=8;
        elseif(bd==32)
            form='single'
            for_mult=4;
        elseif (bd==16)
            form='uint16';
            for_mult=2;
        elseif (bd==8)
            form='uint8';
            for_mult=1;
        end
        for_mult = uint64(for_mult); 
        framenum=num2read;
        imData=cell(1,framenum);
        
        
        % Use low-level File I/O to read the file
        fp = fopen(path_to_file , 'rb');
        % The StripOffsets field provides the offset to the first strip. Based on
        % the INFO for this file, each image consists of 1 strip.
        he_w = info.Width; 
        sizx=he_w;
        he_w=uint64(he_w);
        he_h = info.Height;
        sizy=he_h;
        he_h=uint64(he_h);
        
        first_offset=info.StripOffsets;
        first_offset = uint64(first_offset); 
        ofds=zeros(numFrames, 1, 'like', uint64(0));
        %compute frame offsets
        for i=1:numFrames
            ofds(i)=first_offset+uint64(i-1)*he_w*he_h*for_mult;
            %ofds(i)
        end
        
        
        sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(num2read+sframe-1),' of ',num2str(num_tot_frames), ' total frames'];
        disp(sframemsg)
        pause(.2)
        %go to start of first strip
        %fseek(fp, ofds(1), 'bof');
        % mul is set to > 1 for debugging only
        mul=1;
        
        if strcmpi(form,'uint16') || strcmpi(form,'uint8')
            if(bo)
                for cnt = sframe:lastframe
                    %cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be')';
                    imData{cnt-sframe+1}=cast(tmp1,form);
                end
            else
                for cnt=sframe:lastframe
                    % cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le')';
                    imData{cnt-sframe+1}=cast(tmp1,form);
                end
            end
        elseif strcmpi(form,'single')
            if(bo)
                for cnt = sframe:lastframe
                    %cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be')';
                    imData{cnt-sframe+1}=cast(tmp1,'single');
                end
            else
                for cnt = sframe:lastframe
                    %cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le')';
                    imData{cnt-sframe+1}=cast(tmp1,'single');
                end
            end
        elseif strcmpi(form,'double')
            if(bo)
                for cnt = sframe:lastframe
                    %cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be.l64')';
                    imData{cnt-sframe+1}=cast(tmp1,'single');
                end
            else
                for cnt = sframe:lastframe
                    %cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le.l64')';
                    imData{cnt-sframe+1}=cast(tmp1,'single');
                end
            end
        end
        
        
    else
        blah=size(info);
        numFrames= blah(1);
        num_tot_frames=numFrames;
        
        
        %pretty inelegent kludge to handle imagej big tiff outs (need to pass header info manually in code)
        if (nargin <=3)
            
            %should add more error checking for args... very ugly code below.  works
            %for me after midnight though...
            if nargin<3
                if nargin<2
                    sframe = 1;
                else
                    sframe=cell2mat(varargin(2));
                end
                num2read=numFrames-sframe+1;
            end
            if nargin==3
                num2read=cell2mat(varargin(3));
                sframe = cell2mat(varargin(2));
            end
            if sframe<=0
                sframe=1;
            end
            if num2read<1
                num2read=1;
            end
            
            
            
            if sframe>num_tot_frames
                sframe=num_tot_frames;
                num2read=1;
                display('starting frame has to be less than number of total frames...');
            end
            if (num2read+sframe<= num_tot_frames+1)
                lastframe=num2read+sframe-1;
            else
                num2read=numFrames-sframe+1;
                lastframe=num_tot_frames;
                display('Hmmm...just reading from starting frame until the end');
            end
            
            
            bd=info.BitDepth;
            he=info.ByteOrder;
            bo=strcmp(he,'big-endian');
            if (bd==64)
                form='double';
            elseif(bd==32)
                form='single';
            elseif (bd==16)
                form='uint16';
            elseif (bd==8)
                form='uint8';
            end
            
            
            % Use low-level File I/O to read the file
            fp = fopen(path_to_file , 'rb');
            % The StripOffsets field provides the offset to the first strip. Based on
            % the INFO for this file, each image consists of 1 strip.
            
            he=info.StripOffsets;
            %finds the offset of each strip in the movie.  Image does not have to have
            %uniform strips, but needs uniform bytes per strip/row.
            idss=max(size(info(1).StripOffsets));
            ofds=zeros(numFrames);
            for i=1:numFrames
                ofds(i)=info(i).StripOffsets(1);
                %ofds(i)
            end
            sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(num2read+sframe-1),' of ',num2str(num_tot_frames), ' total frames'];
            disp(sframemsg)
            pause(.2)
            %go to start of first strip
            fseek(fp, ofds(1), 'bof');
            %framenum=numFrames;
            framenum=num2read;
            imData=cell(1,framenum);
            
            he_w=info.Width;
            sizx=he_w;
            he_h=info.Height;
            sizy=he_h;
            % mul is set to > 1 for debugging only
            mul=1;
            if strcmpi(form,'uint16') || strcmpi(form,'uint8')
                if(bo)
                    for cnt = sframe:lastframe
                        %cnt;
                        fseek(fp,ofds(cnt),'bof');
                        tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be')';
                        imData{cnt-sframe+1}=cast(tmp1,form);
                    end
                else
                    for cnt=sframe:lastframe
                        % cnt;
                        fseek(fp,ofds(cnt),'bof');
                        tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le')';
                        imData{cnt-sframe+1}=cast(tmp1,form);
                    end
                end
            elseif strcmpi(form,'single')
                if(bo)
                    for cnt = sframe:lastframe
                        %cnt;
                        fseek(fp,ofds(cnt),'bof');
                        tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be')';
                        imData{cnt-sframe+1}=cast(tmp1,'single');
                    end
                else
                    for cnt = sframe:lastframe
                        %cnt;
                        fseek(fp,ofds(cnt),'bof');
                        tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le')';
                        imData{cnt-sframe+1}=cast(tmp1,'single');
                    end
                end
            elseif strcmpi(form,'double')
                if(bo)
                    for cnt = sframe:lastframe
                        %cnt;
                        fseek(fp,ofds(cnt),'bof');
                        tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be.l64')';
                        imData{cnt-sframe+1}=cast(tmp1,'single');
                    end
                else
                    for cnt = sframe:lastframe
                        %cnt;
                        fseek(fp,ofds(cnt),'bof');
                        tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le.l64')';
                        imData{cnt-sframe+1}=cast(tmp1,'single');
                    end
                end
            end
            
        end
    end
    %ieee-le.l64
    
    imData=cell2mat(imData);
    imData=reshape(imData,[he_h*mul,he_w,framenum]);
    fclose(fp);
    display('Finished reading images')
    
    
elseif strcmpi(ext,'.hdf5') || strcmpi(ext,'.h5');
    info = hdf5info(path_to_file);
    dims = info.GroupHierarchy.Datasets.Dims;
    if nargin < 2
        sframe = 1;
    end
    if nargin < 3
        num2read = dims(end)-sframe+1;
    else
        num2read = round(varargin{3});
    end
    num2read = min(num2read,dims(end)-sframe+1);
    sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(num2read+sframe-1),' of ',num2str(dims(end)), ' total frames'];
    
    imData =squeeze( h5read(path_to_file,info.GroupHierarchy.Datasets.Name,[ones(1,length(dims)-1),sframe],[dims(1:end-1),num2read]));
    display('Finished reading images')
elseif strcmpi(ext,'.avi')
    obj = audiovideo.mmreader(path_to_file);
    
    frame_rate = obj.FrameRate;
    total = obj.Duration;
    numFrames = total*frame_rate;
    
    sizx = obj.Width;
    sizy = obj.Height;
    
    if nargin < 3
        num2read = numFrames-sframe+1;
    else
        num2read = min(numFrames-sframe+1, round(varargin{3}));
    end
    imData = obj.read([sframe, sframe+num2read-1]);
    imData = squeeze(imData(:, :, 1, :));
elseif strcmpi(ext, '.mat')
    data = matfile(path_to_file);
    data_info = whos(data);
    
    if length(data_info)>1
        dims = data.Ysiz;
    else
        dims = data_info.size;
    end
    sizy = dims(1);
    sizx = dims(2);
    numFrames = dims(3);
    if nargin < 3
        num2read = numFrames-sframe+1;
    else
        num2read = min(numFrames-sframe+1, round(varargin{3}));
    end
    
    if length(data_info)>1
        imData = data.Y(:, :, sframe+(0:(num2read-1)));
    else
        imData = eval(sprintf('data.%s(:, :, sframe+(0:(num2read-1)))', data_info.name));
    end
elseif isempty(ext)
    % the input is a folder and data are stored as image sequences
    
    imgs = dir([path_to_file, filesep,'*.tif']);
    frames = zeros(1,length(imgs));
    for file = 1:length(imgs)
        fnam = imgs(file).name;
        IDs = regexp(fnam,'\d*','match'); % use a regular expression to find frame index within each file name
        frames(file) = str2num(IDs{1});
    end
    
    % sort frames in ascending order before reading
    [~,srt] = sort(frames,'ascend');
    imgs = imgs(srt);
    
    %get image info from first image of sequence
    first_img= fullfile(path_to_file,imgs(1).name);
    
    info = imfinfo(first_img);
    sizx = info.Height;
    sizy = info.Width;
    
    num_tot_frames=length(imgs);
    if nargin < 2
        sframe = 1;
    end
    if nargin < 3
        num2read = num_tot_frames-sframe+1;
    else
        num2read = min(num_tot_frames-sframe+1, round(varargin{3}));
    end
    
    if (num2read+sframe<= num_tot_frames+1)
        lastframe=num2read+sframe-1;
    else
        num2read=num_tot_frames-sframe+1;
        lastframe=num_tot_frames;
        display('Hmmm...just reading from starting frame until the end');
    end
    
    imData = zeros(sizx,sizx,num2read);
    
    sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(lastframe),' of ',num2str(num_tot_frames), ' total frames'];
    disp(sframemsg)
    
    for t = 1:num2read
        imData(:,:,t)=imread(fullfile(path_to_file,imgs(t+sframe-1).name));
    end
    
    display('Finished reading images')
else
    error('Unknown file extension. Only .tiff and .hdf5 files are currently supported');
end
