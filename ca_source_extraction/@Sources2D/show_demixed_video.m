function show_demixed_video(obj,save_avi, kt, center_ac, range_ac, range_Y, multi_factor)
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
%% load data 
frame_range = obj.frame_range; 
t_begin = frame_range(1); 
t_end = frame_range(2); 
%% data preparation
Y = obj.load_patch_data(); 
Ybg = obj.reconstruct_background();
figure('position', [0,0, 600, 400]);

%%
if ~exist('center_ac', 'var') || isempty(center_ac)
    center_ac = median(max(obj.A,[],1)'.*max(obj.C,[],2));
end
range_res = [-1,1]*center_ac;
if ~exist('range_ac', 'var') || isempty(range_ac)
    range_ac = center_ac*1.01+range_res;
end
if ~exist('range_Y', 'var') || isempty(range_Y)
    if ~exist('multi_factor', 'var') || isempty(multiple_factor)
        temp = quantile(double(Y(randi(numel(Y), 10000,1))), [0.01, 0.98]);
        multi_factor = ceil(diff(temp)/diff(range_ac));
    else
        temp = quantile(Y(randi(numel(Y), 10000,1)), 0.01); 
    end
    center_Y = temp(1) + multi_factor*center_ac;
    range_Y = center_Y + range_res*multi_factor;
end
%% create avi file
if save_avi
    avi_filename =[obj.P.log_folder, 'demixed.avi'];
    avi_file = VideoWriter(avi_filename);
    if ~isnan(obj.Fs)
        avi_file.FrameRate= obj.Fs/kt;
    end
    avi_file.open();
end

%% add pseudo color to denoised signals
[K, T]=size(obj.C);
% draw random color for each neuron
% tmp = mod((1:K)', 6)+1;
Y_mixed = zeros(obj.options.d1*obj.options.d2, T, 3);
temp = prism;
% temp = bsxfun(@times, temp, 1./sum(temp,2));
col = temp(randi(64, K,1), :);
for m=1:3
    Y_mixed(:, :, m) = obj.A* (diag(col(:,m))*obj.C);
end
Y_mixed = uint16(Y_mixed/(1*center_ac)*65536);
%% play and save
ax_y =   axes('position', [0.015, 0.51, 0.3, 0.42]);
ax_bg=   axes('position', [0.015, 0.01, 0.3, 0.42]);
ax_signal=    axes('position', [0.345, 0.51, 0.3, 0.42]);
ax_denoised =    axes('position', [0.345, 0.01, 0.3, 0.42]);
ax_res =    axes('position', [0.675, 0.51, 0.3, 0.42]);
ax_mix =     axes('position', [0.675, 0.01, 0.3, 0.42]);
for m=t_begin:kt:t_end
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
    img_ac = obj.reshape(obj.A*obj.C(:, m), 2); 
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
    text(1, 10, sprintf('Time: %.2f second', m/obj.Fs), 'color', 'w', 'fontweight', 'bold');
    
    axis equal tight off;
    %     box on; set(gca, 'xtick', []);
    %     set(gca, 'ytick', []);
    
    drawnow(); 
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [400, 600]);
        avi_file.writeVideo(temp);
    end
end

if save_avi
    avi_file.close();
end
