%%
Yac = neuron.reshape(neuron.A*neuron.C, 2);
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
    avi_file = VideoWriter([dir_nm, filesep, file_nm, '_results.avi']);
<<<<<<< Updated upstream
    avi_file.FrameRate = neuron.Fs; 
=======
    avi_file.FrameRate= neuron.Fs/kt; 
>>>>>>> Stashed changes
    avi_file.open();
end
temp  = quantile(Ybg(1:1000:numel(Ybg)), [0.0001, y_quantile]);
Ymin = temp(1);
Ymax = temp(2);
Ymin_1 = temp(1) - mean(Ymean(:));
Ymax_1 = temp(2) - mean(Ymean(:));

ACmax = quantile(Yac(1:1000:numel(Yac)), ac_quantile);
ACmin = ACmax/100; 

%     subplot(4,6, [5,6,11,12]);
<<<<<<< Updated upstream
for m=5381:kt:T
=======
for m=t_begin:kt:t_end
>>>>>>> Stashed changes
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
    
    
    subplot(4, 6, [3,4, 9,10]+12);
    imagesc(Ysignal(:, :, m), [ACmin, ACmax]); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw-background'); hold on; colorbar;
    
        
    subplot(4,6, [5, 6, 11, 12]+12);    
    imagesc(Yac(:, :, m), [ACmin, ACmax]);
    %     imagesc(Ybg(:, :, m), [-50, 50]);
    text(d2/2-30, 10, sprintf('Time: %.2f Sec', m/neuron.Fs), 'color', 'w');
    axis equal off tight; title('denoised'); colorbar;
    
subplot(4,6, [5,6,11,12]);
    imagesc(Ysignal(:, :, m)-Yac(:, :, m), [-ACmax, ACmax]/2); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
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
