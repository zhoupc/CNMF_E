%% change gaussian width 
tmp_gsig = str2double(get(edit_gsig, 'string')); 
neuron_raw.options.gSig = tmp_gsig; 
if ~exist('neuron', 'var')
    neuron.options.gsig = tmp_gsig / neuron_raw.options.ssub; 
end