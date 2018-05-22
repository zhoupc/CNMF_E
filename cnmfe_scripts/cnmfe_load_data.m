if and(ssub==1, tsub==1)
    neuron = neuron_full;

    tmpY = data.Y; 

    Y = single(tmpY(:, :, sframe+(1:num2read)-1));
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else
    if exist(nam_mat, 'file')
        [Y, neuron_ds] = neuron_full.load_data(nam_mat, sframe, num2read);
        [d1s,d2s, T] = size(Y);
        fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
        neuron = neuron_ds.copy();
    else
        [Y, neuron_ds] = neuron_full.load_data(nam, sframe, num2read);
        [d1s,d2s, T] = size(Y);
        fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
        neuron = neuron_ds.copy();
    end
end
