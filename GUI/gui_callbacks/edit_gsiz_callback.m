%% change gaussian width 
tmp_gsiz = str2double(get(edit_gsiz, 'string')); 
neuron_raw.options.gSiz = tmp_gsiz; 
if ~exist('neuron', 'var')
    neuron.options.gsiz = tmp_gsiz / neuron_raw.options.ssub; 
end