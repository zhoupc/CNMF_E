tic;
neuron.options.nk = 1; %round(T/(60*neuron.Fs)); % number of knots for spline basis, the interval between knots is 180 seconds
nr_patch = str2double(get(edit_nr, 'string')); 
nc_patch = str2double(get(edit_nc, 'string')); 

patch_par = round([nr_patch, nc_patch]); %1;  % divide the optical field into 3 X 3 patches and do initialization patch by patch
K = str2double(get(edit_K, 'string'));  % maximum number of neurons to search within each patch. you can use [] to search the number automatically
neuron.options.bd = 0; % boundaries to be removed due to motion correction
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); %#ok<ASGLU>
figure;
imagesc(Cn);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;

[~, srt] = sort(max(neuron.C, [], 2)./get_noise_fft(neuron.C), 'descend');
neuron.orderROIs(srt);
neuron_init = neuron.copy();
