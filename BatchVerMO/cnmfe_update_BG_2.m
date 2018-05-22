clear Ysignal;
tic;
display(size(NEURON.A))
display(size(NEURON.C))
display(size(Y))
Ybg = Y-NEURON.A*NEURON.C;
rr = ceil(NEURON.options.gSiz * bg_neuron_ratio); 
active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
[Ybg, Ybg_weights] = NEURON.localBG(Ybg, spatial_ds_factor, rr, active_px, NEURON.P.sn, thresh); % estiamte local background.
% subtract the background from the raw data.
Ysignal = Y - Ybg;

%% estimate noise 
if ~isfield(NEURON.P, 'sn') || isempty(NEURON.P.sn)
    %% estimate the noise for all pixels
    b0 =zeros(size(Ysignal,1), 1);
    sn = b0;
    parfor m=1:size(NEURON.A,1)
        [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m, :));
    end     
    Ysignal = bsxfun(@minus, Ysignal, b0);
    NEURON.P.sn = sn; 
end