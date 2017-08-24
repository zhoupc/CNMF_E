function img = spatial_constraints(img, show_imgs)
% enforce a stronger prior to neuron's spatial component by assuming that
% the gradient at each pixel is pointing to the peak location 
% It also removes isolated pixels. 
if nargin<2
    show_imgs = false; 
end
[tmp1, tmp2, ~] = find(img); 
rmin = min(tmp1); 
rmax = max(tmp1); 
cmin = min(tmp2); 
cmax = max(tmp2); 
[nr, nc] = size(img);   % image dimension  

% crop a small region
if rmin==1 && rmax==nr && cmin==1 && cmax==nc
    
    if show_imgs
        figure;
        subplot(121); 
        imagesc(img);
        axis equal off tight;
    end
    
    [vmax, ind_max] = max(img(:));
    [y0, x0] = ind2sub([nr, nc], ind_max);
    [x, y] = meshgrid(1:nc, 1:nr);
    [fx, fy] = gradient(img);
    ind = ((fx.*(x0-x)+fy.*(y0-y)) < 0);% & (img<vmax/2);
    img(ind) = 0;
    l = bwlabel(img, 4);
    img(l~=l(ind_max)) = 0;
    img = medfilt2(img); 
    
    if show_imgs
        subplot(122);
        imagesc(img);
        axis equal off tight;
        pause; 
        close; 
    end
else
    tmp_img = spatial_constraints(img(rmin:rmax, cmin:cmax), show_imgs);
    img(rmin:rmax, cmin:cmax) = tmp_img;
end

