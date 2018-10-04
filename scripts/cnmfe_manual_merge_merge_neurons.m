% determine searching area
active_pixel = (sum(neuron_bk.A(:,IDs), 2)>0);

% update spatial/temporal components of the merged neuron
data = neuron_bk.A(active_pixel, IDs)*neuron_bk.C_raw(IDs, :);
ci = neuron_bk.C_raw(IDs(1), :);
for miter=1:10
    ai = data*ci'/(ci*ci');
    ci = ai'*data/(ai'*ai);
end

sn = GetSn(ci);
neuron.A(active_pixel, IDs(1)) = ai*sn;
neuron.C_raw(IDs(1), :) = ci/sn;
neuron.S(IDs(1), :) = ci/sn; 
newIDs(IDs(1)) = IDs(1);
% remove merged elements
ind_del(IDs(2:end)) = true;

%% show merged shape and trace 
axes(ax_merged); 
neuron.image(neuron.A(:, IDs(1))); 
axis equal off tight; 
set(gca, 'xlim', get(ax_selected, 'xlim')); 
set(gca, 'ylim', get(ax_selected, 'ylim')); 

axes(ax_merged_trace); 
plot(ci, 'k', 'linewidth', 2); 

%% clear IDs; 
IDs_last = IDs; 
IDs = []; 