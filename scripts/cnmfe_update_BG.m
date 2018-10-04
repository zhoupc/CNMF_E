clear Ysignal;
tic; 
Ybg = Y-neuron.A*neuron.C;
rr = ceil(neuron.options.gSiz * bg_neuron_ratio); 
active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
[Ybg, Ybg_weights] = neuron.localBG(Ybg, spatial_ds_factor, rr, active_px, neuron.P.sn, thresh); % estiamte local background.
% subtract the background from the raw data.
Ysignal = Y - Ybg;

%% estimate noise 
if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
    %% estimate the noise for all pixels
    b0 =zeros(size(Ysignal,1), 1);
    sn = b0;
    parfor m=1:size(neuron.A,1)
        [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m, :));
    end     
    Ysignal = bsxfun(@minus, Ysignal, b0);
    neuron.P.sn = sn; 
end