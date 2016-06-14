%% callback of after selecting searching method
if get(radio_dilate, 'value');
    set(radio_ellipse, 'value', 0);
    bSiz = str2double(get(edit_dilate, 'string'));
    
    if exist('neuron', 'var')
        neuron.options.search_method = 'dilate';
        neuron.options.se = strel('disk', bSiz);
        neuron.options.bSiz = bSiz;
    else
        neuron_raw.options.search_method = 'dilate';
        neuron_raw.options.se = strel('disk', bSiz);
        neuron_raw.options.bSiz = bSiz;
    end
else
    set(radio_ellipse, 'value', 1);
    radio_ellipse_callback;
end