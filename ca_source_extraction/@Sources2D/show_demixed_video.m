function avi_filename = show_demixed_video(obj,save_avi, kt, frame_range, amp_ac, range_ac, range_Y, multi_factor, use_craw)
%% save the final results of source extraction.
%% inputs:
%   kt:  scalar, the number of frames to be skipped

%% author: Pengcheng Zhou, Columbia University 2017
% zhoupc1988@gmail.com
if ~exist('save_avi', 'var')||isempty(save_avi)
    save_avi = false;
end
if ~exist('kt', 'var')||isempty(kt)
    kt = 1;
end
if ~exist('use_craw', 'var')||isempty(use_craw)
    use_craw = false; 
end
if use_craw
    CC = obj.C_raw; 
else
    CC = obj.C; 
end

%% load data
if ~exist('frame_range', 'var')||isempty(frame_range)
    frame_range = obj.frame_range;
end
t_begin = frame_range(1);
t_end = frame_range(2);

tmp_range = [t_begin, min(t_begin+100*kt-1, t_end)];
Y = obj.load_patch_data([], tmp_range);
Ybg = obj.reconstruct_background(tmp_range);
figure('position', [0,0, 600, 400]);

%%
if ~exist('amp_ac', 'var') || isempty(amp_ac)
    amp_ac = median(max(obj.A,[],1)'.*max(CC,[],2))*2;
end
if ~exist('range_ac', 'var') || isempty(range_ac)
    range_ac = amp_ac*[0.01, 1.01];
end
range_res = range_ac - mean(range_ac); 
if ~exist('range_Y', 'var') || isempty(range_Y)
    if ~exist('multi_factor', 'var') || isempty(multi_factor)
        temp = quantile(double(Y(randi(numel(Y), 10000,1))), [0.01, 0.98]);
        multi_factor = ceil(diff(temp)/diff(range_ac));
    else
        temp = quantile(Y(randi(numel(Y), 10000,1)), 0.01);
    end
    center_Y = temp(1) + multi_factor*amp_ac/2;
    range_Y = center_Y + range_res*multi_factor;
end

if ~exist('multi_factor', 'var')
    multi_factor = round(diff(range_Y)/diff(range_ac)); 
    range_Y = (range_res-range_res(1)) * multi_factor +range_Y(1); 
end
%% create avi file
if save_avi
    avi_filename =[obj.P.log_folder, 'demixed.avi'];
    avi_file = VideoWriter(avi_filename);
    if ~isnan(obj.Fs)
        avi_file.FrameRate= obj.Fs/kt;
    end
    avi_file.open();
else
    avi_filename = [];
end

%% add pseudo color to denoised signals
[K, ~]=size(CC);
% draw random color for each neuron
% tmp = mod((1:K)', 6)+1;
Y_mixed = zeros(obj.options.d1*obj.options.d2, diff(tmp_range)+1, 3);
temp = prism;
% temp = bsxfun(@times, temp, 1./sum(temp,2));
col = temp(randi(64, K,1), :);
for m=1:3
    Y_mixed(:, :, m) = obj.A* (diag(col(:,m))*CC(:, tmp_range(1):tmp_range(2)));
end
Y_mixed = uint16(Y_mixed*2/(amp_ac)*65536);
%% play and save
ax_y =   axes('position', [0.015, 0.51, 0.3, 0.42]);
ax_bg=   axes('position', [0.015, 0.01, 0.3, 0.42]);
ax_signal=    axes('position', [0.345, 0.51, 0.3, 0.42]);
ax_denoised =    axes('position', [0.345, 0.01, 0.3, 0.42]);
ax_res =    axes('position', [0.675, 0.51, 0.3, 0.42]);
ax_mix =     axes('position', [0.675, 0.01, 0.3, 0.42]);
tt0 = (t_begin-1);
for tt=t_begin:kt:t_end
    m = tt-tt0;
    
    axes(ax_y); cla;
    imagesc(Y(:, :,m), range_Y);
    
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title('Raw data');
    axis equal off tight;
    
    axes(ax_bg); cla;
    imagesc(Ybg(:, :, m),range_Y);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    axis equal off tight;
    title('Background');
    
    axes(ax_signal); cla;
    imagesc(double(Y(:, :, m))-Ybg(:, :, m), range_ac); hold on;
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title(sprintf('(Raw-BG) X %d', multi_factor));
    axis equal off tight;
    
    axes(ax_denoised); cla;
    img_ac = obj.reshape(obj.A*CC(:, tt), 2);
    imagesc(img_ac, range_ac);
    %     imagesc(Ybg(:, :, m), [-50, 50]);
    title(sprintf('Denoised X %d', multi_factor));
    axis equal off tight;
    
    axes(ax_res); cla;
    imagesc(double(Y(:, :, m))-Ybg(:, :, m)-img_ac, range_res);
    %     set(gca, 'children', flipud(get(gca, 'children')));
    title(sprintf('Residual X %d', multi_factor));
    axis equal off tight;
    %         subplot(4,6, [5,6,11,12]+12);
    
    axes(ax_mix); cla;
    imagesc(obj.reshape(Y_mixed(:, m,:),2));  hold on;
    title('Demixed');
    text(1, 10, sprintf('Time: %.2f second', (tt)/obj.Fs), 'color', 'w', 'fontweight', 'bold');
    
    axis equal tight off;
    %     box on; set(gca, 'xtick', []);
    %     set(gca, 'ytick', []);
    
    drawnow();
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [400, 600]);
        avi_file.writeVideo(temp);
    end
    
    %% save more data
    if tt== tt0+kt*99+1
        tt0 = tt0+kt*100;
        % load data
        tmp_range = [tt0+1, min(tt0+100*kt, t_end)];
        if isempty(tmp_range) || diff(tmp_range)<1
            break;
        end
        Y = obj.load_patch_data([], tmp_range);
        
        Ybg = obj.reconstruct_background(tmp_range);
        
        %% add pseudo color to denoised signals
        [d1, d2, Tp]=size(Y);
        % draw random color for each neuron
        % tmp = mod((1:K)', 6)+1;
        Y_mixed = zeros(d1*d2, Tp, 3);
        tmp_C = CC(:, tt0+(1:Tp));
        temp = prism;
        col = temp(randi(64, K,1), :);
        for m=1:3
            Y_mixed(:, :, m) = obj.A* (diag(col(:,m))*tmp_C);
        end
        Y_mixed = uint16(Y_mixed*2/(amp_ac)*65536);
    end
end

if save_avi
    avi_file.close();
end
