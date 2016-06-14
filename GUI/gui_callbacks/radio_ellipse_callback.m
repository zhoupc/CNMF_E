%% callback of after selecting searching method
if get(radio_ellipse, 'value')
    set(radio_dilate, 'value', 0);
    dist_ratio = str2double(get(edit_ellipse, 'string'));
    
    if exist('neuron', 'var')
        neuron.options.search_method = 'ellipse';
        neuron.options.dist = dist_ratio;
    else
        neuron_raw.options.search_method = 'ellipse';
        neuron_raw.options.dist = dist_ratio;
    end
else
    set(radio_dilate, 'value', 1);
    radio_dilate_callback;
end
