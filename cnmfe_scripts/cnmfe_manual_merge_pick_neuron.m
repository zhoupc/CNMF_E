%% manually pick neurons
axes(ax_all);
xrange = get(gca, 'xlim');
yrange = get(gca, 'ylim');
title('click pixels out of the field to end');
while true
    [tmp_x, tmp_y] = ginput(1);
    if tmp_x<xrange(1) || tmp_x>xrange(2) || (tmp_y<yrange(1)) || (tmp_y>yrange(2))
        break;
    else
        % find the closest neuron
        dist = (tmp_x- ctr(:,2)).^2 + (tmp_y - ctr(:,1)).^2;
        [~, tmp_ID] = min(dist);
        tmp_col = color_all(tmp_ID, :);
        plot(tmp_x, tmp_y, '*', 'color', tmp_col);
        if any(IDs==tmp_ID) || ind_del(tmp_ID) % already been selected or merged
            continue;
        else
            IDs = [IDs; tmp_ID];
            plot(Coor{tmp_ID}(1,:), Coor{tmp_ID}(2,:), 'color', tmp_col);
        end
    end
end

% compute distance to all neurons
if isempty(IDs)
    return;
end

%% show spatial shapes of the selected neurons
img_selected = zeros(size(AA,1), 3);
col = col0;
for m=1:3
    img_selected(:, m) = sum(bsxfun(@times, AA(:, IDs), mod(col(:, IDs), 2)), 2);
    col = floor(col/2);
end
img_selected = neuron.reshape(img_selected, 2);
img_selected = img_selected/max(img_selected(:))*(2^16);
img_selected = uint16(img_selected);

axes(ax_selected);
imagesc(img_selected);
axis equal off tight;
xmin = min(ctr(IDs, 2));
xmax = max(ctr(IDs, 2));
ymin = min(ctr(IDs, 1));
ymax = max(ctr(IDs, 1));
gSiz = neuron.options.gSiz;
xlim([xmin-gSiz, xmax+gSiz]);
ylim([ymin-gSiz, ymax+gSiz]);

%% show temporal traces of the selected neurons
tmp_C = neuron.C_raw(IDs, :)';
tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
axes(ax_trace);cla;  hold on;
for m=1:length(IDs)
    plot(tmp_C(:,m), 'color', color_all(IDs(m),:),  'linewidth', 2);
end







