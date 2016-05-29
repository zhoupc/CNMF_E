clear; clc;
%% choose data
if ~exist('nam', 'var') || isempty(nam)
    try
        load .dir.mat;
    catch
        dir_nm = [cd(), filesep];
    end
    [file_nm, dir_nm] = uigetfile(sprintf('%s*.tif', dir_nm));
    save .dir.mat dir_nm;
    nam = [dir_nm, file_nm];
end

temp = imread(nam);
[d1,d2, ~] = size(temp);
neuron_ica_raw = Sources2D('d1',d1,'d2',d2);                         % dimensions of datasets
neuron_ica_raw.updateParams('ssub', 1, 'tsub', 1);    % it seems that the
% raw data has been downsampled by ICA. adjust the downsampling factor
% here

%% load data
tic;
sframe=1;						% user input: first frame to read (optional, default 1)
num2read= 9000;					% user input: how many frames to read   (optional, default until the end)
[Y, neuron_ica] = neuron_ica_raw.load_data(nam, sframe, num2read);
[d1,d2, T] = size(Y);
Y = neuron_ica.reshape(Y, 1);
fprintf('Time cost in loading data:     %.2f seconds\n', toc);

%% directory to pca/ica results
dir_ics = uigetdir(dir_nm, 'folder of all ICs');
dir_traces = uigetdir(fileparts(dir_ics), 'folder of the traces');
dir_events = uigetdir(fileparts(dir_traces), 'folder of the events');

%% load calcium traces
C = load(sprintf('%s%straces.txt', dir_traces, filesep));

t = C(1:T, 1);
C = C(1:T, 2:end)';
neuron_ica.C = C;
K = size(C, 1);
% load events
S = zeros(K, T);
temp = load(sprintf('%s%sevents.txt', dir_events, filesep));
S(1, :) = temp(1:T, 2);
for m=2:K
    temp = load(sprintf('%s%sevents%d.txt', dir_events, filesep, m-1));
    S(m, :) = temp(1:T, 2);
end
neuron_ica.S = S;
% load spatial components
A = zeros(d1, d2, K);
ssub = neuron_ica.options.ssub;
A(:, :, 1) = imresize(imread(sprintf('%s%sics.tif', dir_ics,filesep)), [d1,d2]);
for m=2:K
    A(:, :, m) = imresize(imread(sprintf('%s%sics%d.tif', dir_ics, filesep, m-1)), [d1,d2]);
end
neuron_ica.A = neuron_ica.reshape(A, 1);


%%
% min_max = [500, 2000]; % image scaling
% save_avi = false;
% avi_name = 'ica_check.avi';
% neuron_ica.runMovie(Y, min_max, save_avi, avi_name, S);

%% compute a data for fast initialization
psf = fspecial('gaussian', 13, 3); % this gaussian kernel is a rough estimation of neuron shape 
% psf = psf-mean(psf(:));
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;
HY = imfilter(neuron_ica.reshape(Y, 2), psf, 'replicate');
HY = neuron_ica.reshape(HY, 1);
% compute local correlatin image 
dHY = diff(HY(:, 1:3:end), 1, 2); 
dHY(bsxfun(@lt, dHY, 2*std(dHY, 0, 2))) = 0; 
Cn = correlation_image(dHY, 8, d1, d2); 
neuron.Cn = Cn; 

%% plot contours
neuron_ica.viewContours(Cn, 1, 0);

%%
figure('position', [100, 100, 1800, 420]);

% display correlation image for selecting neurons
subplot(2, 5, [1,2,6,7]);
clim = [0.3, 1]; 
neuron_ica.image(Cn, clim); hold on;
axis equal; axis off tight;hold on;
m = 1;
ctr = neuron_ica.estCenter();
Coor = neuron_ica.Coor;
for m=1:length(Coor)
    % remove those weird shape contours
    temp = Coor{m};
    ind = or(abs(temp(1,:)-ctr(m,2))>10, abs(temp(2,:)-ctr(m,1))>10);
    temp(:, ind) = [];
    Coor{m} = temp;
    plot(temp(1, :), temp(2, :), 'r');
end
gSiz = 20;  % window size of the cropped video 
Fr = 100;  % frame rate when playing the video 
tstart = 1;  
tend = T;
kk=1;
while true
    subplot(2,5,[1,2,6,7]);
    try  %#ok<TRYNC>
        delete(h_con);
        delete(h_txt);
        delete(h_dot);
    end
    [x, y] = ginput(1);
%     kk = kk+1;
%     y = ctr(kk, 1);
%     x = ctr(kk, 2);
    x = round(x); y = round(y);
    h_dot = plot(x, y, '*m');
    dx = x-ctr(:, 2);
    dy = y-ctr(:, 1);
    dist_xy = sqrt(dx.^2+dy.^2);
    [~, ind_cell] = min(dist_xy);
    temp = Coor{ind_cell};
    h_con = plot(temp(1, :), temp(2, :), 'g', 'linewidth', 2);
    h_txt = text(ctr(ind_cell, 2), ctr(ind_cell, 1), num2str(ind_cell), 'fontsize', 10);
    
    axis equal off tight;
    subplot(2,5,3:5); cla;
    plot(C(ind_cell, :)); hold on;
    
    if (exist('S', 'var'))&&  (~isempty(S))
        tmp = find(S(ind_cell, :)>0);
        plot(tmp, S(ind_cell, tmp), '*g');
    end
    set(gca, 'xticklabel', []); axis tight;
    xlim([1, T]);
    title('ICA result');
    subplot(2,5,8:10); cla;
    plot(squeeze(HY(sub2ind([d1,d2], round(y), round(x)), :))); hold on;
    title('trace on the purple pixel');
    xlim([1, T]);
    set(gca, 'xticklabel', get(gca, 'xtick'));
    xlabel('Frame');
    
    choice = questdlg('what'' next?', 'make a decision', 'continue', ...
        'play video', 'stop', 'continue');
    switch choice
        case 'play video'
            % input some parameters
            prompt = {'window size', 'time start', 'time end', 'frame rate', 'save video'};
            defautans = {num2str(gSiz), num2str(tstart), num2str(tend), num2str(Fr), 'File Name'};
            answer = inputdlg(prompt, 'parameters for playing video', [1, 20], defautans);
            gSiz = round(str2double(answer{1}));
            tstart = round(str2double(answer{2}));
            tend = round(str2double(answer{3}));
            Fr = round(str2double(answer{4}));
            file_nm = answer{5};
            if ~strcmp('File Name', file_nm)
                fp = VideoWriter([file_nm, '.avi']); %#ok<TNMLP>
                fp.FrameRate = Fr;
                fp.open();
                frame = getframe(gcf);
                for ff=1:Fr
                    fp.writeVideo(frame);
                end
            end
            
            [indc, indr] = meshgrid(max(1,x-gSiz):min(d2, x+gSiz), max(1, y-gSiz):min(d1, y+gSiz));
            [nr, nc] = size(indc);
            r0 = y-indr(1);
            c0 = x-indc(1);
            tmp_cont = bsxfun(@minus, Coor{ind_cell}, [indc(1), indr(1)]');
            % crop the video 
            data = Y(sub2ind([d1,d2], indr, indc), :);
            % estimate background using boundary pixels 
            indr = [ones(1, nc), ones(1, nc)*nr, 1:nr, 1:nr];
            indc = [1:nc, 1:nc, ones(1, nr)*nc, ones(1,nr)];
            ind_bd = sub2ind([nr, nc], indr, indc);     % indices for boundary pixels
            y_bg = median(Y(ind_bd, :), 1);  % take the mean in the boundary as an estimation of the background level

            f0 = mean(data(Cn(sub2ind([d1,d2], indr, indc))<3, :), 1);
            data = data - min(bsxfun(@times, data, 1./f0), [], 2)*f0;
            data = reshape(data, nr, nc, []);
            subplot(2,5,3:5); hold on;
            a1 = plot([1,1]*tstart, get(gca, 'ylim'), 'r');
            subplot(2,5,8:10); hold on;
            a2 = plot([1,1]*tstart, get(gca, 'ylim'), 'r');
            subplot(2,5,[1,2,6,7]); cla; hold on;
            tmp_clim = quantile(data(:), [0.05, 1]);
            plot(c0, r0, '*m');
            plot(tmp_cont(1,:), tmp_cont(2,:), 'g', 'linewidth',2);
            for m=tstart:tend
                if m==tstart
                    tic;
                end
                subplot(2,5,[1,2,6,7]);
                if m~=tstart
                    delete(tmp_img);
                end
                tmp_img = imagesc(data(:, :, m), tmp_clim);
                
                set(gca, 'children', flipud(get(gca, 'children')));
                axis equal off tight;
                set(a1, 'XData', [m, m]);
                set(a2, 'Xdata', [m, m]);
                if m==tstart
                    t_plot = toc();
                end
                drawnow();
                pause(max(0.001, 1/Fr-t_plot));
                if ~strcmp('File Name', file_nm)
                    frame = getframe(gcf);
                    fp.writeVideo(frame);
                end
            end
            if ~strcmp('File Name', file_nm)
                fp.close();
            end
            subplot(2,5,[1,2,6,7]); cla;
            neuron_ica.image(Cn, clim);
            hold on;
            for m=1:length(Coor)
                temp = Coor{m};
                plot(temp(1, :), temp(2, :), 'r');
            end
        case 'stop'
            break;
        case 'continue'
            continue;
    end
end















