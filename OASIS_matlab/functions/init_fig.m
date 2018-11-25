set(gcf, 'color', 'w', ...
    'defaultAxesFontSize', 20, ...)
    'defaultlinelinewidth',2, ...
    'paperposition', [0, 0, get(gcf, 'papersize')]); 

% screen size 
tmp_fig_size = round(100*get(gcf, 'papersize')); 
screen_size = get(0, 'ScreenSize'); 
if ~exist('fig_x0', 'var')
    fig_x0 = screen_size(3)/2 - tmp_fig_size(1)/2;
end
if ~exist('fig_y0', 'var')
    fig_y0 = screen_size(4)/2 - tmp_fig_size(2)/2;
end
set(gcf, 'position', [fig_x0, fig_y0, tmp_fig_size]); 

clear fig_x0 fig_y0; 