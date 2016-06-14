%% comp%% compute correlation image and peak-to-noise ratio image.
tmp_h = msgbox('Please wait, I am working...'); 
[Cn, pnr] = neuron.correlation_pnr(Y);
delete(tmp_h); 
figure('position', [10, 500, 1776, 400]);
subplot(131);
imagesc(Cn, [0, 1]); colorbar;
axis equal off tight;
title('correlation image');

subplot(132);
imagesc(pnr,[0,max(pnr(:))*0.98]); colorbar;
axis equal off tight;
title('peak-to-noise ratio');

subplot(133);
imagesc(Cn.*pnr, [0,max(pnr(:))*0.98]); colorbar;
axis equal off tight
title('Cn*PNR'); 