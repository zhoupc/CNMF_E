function [xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size,gridstartend)
if isempty(gridstartend)
    start_xx=1; start_yy=1; start_zz=1;
    end_xx=d1;  end_yy=d2;  end_zz=d3;
else
    start_xx=gridstartend(1); end_xx=gridstartend(2); 
    start_yy=gridstartend(3); end_yy=gridstartend(4);
    start_zz=gridstartend(5); end_zz=gridstartend(6);
end
xx_s = start_xx:grid_size(1):end_xx;
yy_s = start_yy:grid_size(2):end_yy;
zz_s = start_zz:grid_size(3):end_zz;

xx_f = [xx_s(2:end)-1,end_xx];
yy_f = [yy_s(2:end)-1,end_yy];

if length(zz_s)>1
    zz_f = [zz_s(2:end)-1,end_zz];
else
    zz_f=end_zz;
end

if xx_f(end)-xx_s(end) + 1 < min_patch_size(1) && length(xx_s) > 1; xx_s(end) = []; xx_f(end-1) = []; end
if yy_f(end)-yy_s(end) + 1 < min_patch_size(2) && length(yy_s) > 1; yy_s(end) = []; yy_f(end-1) = []; end
if zz_f(end)-zz_s(end) + 1 < min_patch_size(3) && length(zz_s) > 1; zz_s(end) = []; zz_f(end-1) = []; end

grid_size_us = floor(grid_size./mot_uf);
if mot_uf(1) > 1
    xx_us = start_xx:grid_size_us(1):end_xx;
    xx_uf = [xx_us(2:end)-1,end_xx];
else
    xx_us = xx_s; xx_uf = xx_f;
end
if mot_uf(2) > 1
    yy_us = start_yy:grid_size_us(2):end_yy;
    yy_uf = [yy_us(2:end)-1,end_yy];
else
    yy_us = yy_s; yy_uf = yy_f;
end
if mot_uf(3) > 1
    zz_us = start_zz:grid_size_us(3):end_zz;
    zz_uf = [zz_us(2:end)-1,end_zz];
else
    zz_us = zz_s; zz_uf = zz_f;
end