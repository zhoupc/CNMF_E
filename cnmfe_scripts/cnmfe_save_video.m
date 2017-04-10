%%
Yac = neuron.reshape(neuron.A*neuron.C_raw, 2);
Ymean = neuron.reshape(mean(Y, 2), 2);  % estimate the baseline
Ybg = neuron.reshape(Ybg, 2);

% Y = neuron.reshape(Y, 2);
Ysignal = neuron.reshape(Ysignal, 2);
% ctr = round( neuron.estCenter());
% figure;
% neuron.viewContours(Cn, .5, 0);
% cell_IDs = [];

figure('position', [0,0, 1248, 600]);
% avi_file = VideoWriter('~/Dropbox/Public/Garret/day1/residual.avi');
if save_avi
    if ~exist('avi_filename', 'var')
        avi_filename =[dir_nm, filesep, file_nm];
    end
    avi_file = VideoWriter(avi_filename);
    avi_file.FrameRate= neuron.Fs/kt;
    avi_file.open();
end

if ~exist('Ymin', 'var')
    y_quantile = 0.9999;    % for specifying the color value limits
    ac_quantile = .9999;
    temp  = quantile(Ybg(randi(numel(Ybg(:)), 10000,1)), [0.0001, y_quantile]);
    Ymin = temp(1);
    Ymax = temp(2);
    Ymax_1 = Ymax - min(Ymean(:)); 
    ACmax = quantile(Yac(randi(numel(Yac(:)), 10000,1)), ac_quantile);
    ACmin = ACmax/10;
end
%     subplot(4,6, [5,6,11,12]);
for m=t_begin:kt:t_end
    subplot(4,6, [1,2, 7, 8]);
    imagesc(Ybg(:, :,m)+Ysignal(:, :, m), [Ymin, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw data'); hold on; colorbar;
    
    subplot(4, 6, [1,2, 7, 8]+12);
    imagesc(Ybg(:, :, m), [Ymin, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('background');
    colorbar;
    
    subplot(4,6, [3,4, 9,10]);
    imagesc(Ybg(:, :,m)+Ysignal(:, :, m)-Ymean, [0, Ymax_1]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw - mean'); hold on; colorbar;
    
    
    subplot(4,6, [5,6,11,12]);
    imagesc(Ysignal(:, :, m), [ACmin, ACmax]); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw-background'); hold on; colorbar;
    
    
    subplot(4, 6, [3,4, 9,10]+12);
    imagesc(Yac(:, :, m), [ACmin, ACmax]);
    %     imagesc(Ybg(:, :, m), [-50, 50]);
    axis equal off tight; title('denoised'); colorbar;
    
    subplot(4,6, [5,6,11,12]+12);
    imagesc(Ysignal(:, :, m)-Yac(:, :, m), [-ACmax, ACmax]/2); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    text(1, 10, sprintf('Time: %.2f second', m/neuron.Fs), 'color', 'w');
    axis equal; axis off tight; title('residual'); hold on; colorbar;
    
    drawnow();
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [600, 1248]);
        avi_file.writeVideo(temp);
    end
end

if save_avi
    avi_file.close();
end
