neuron.A(:, IDs_last) = neuron_bk.A(:, IDs_last); 
neuron.C(IDs_last,:) = neuron_bk.C(IDs_last,:); 
neuron.C_raw(IDs_last,:) = neuron_bk.C_raw(IDs_last,:); 
ind_del(IDs_last) = false; 

axes(ax_merged); cla; 
axes(ax_merged_trace); cla; 

IDs = IDs_last; 