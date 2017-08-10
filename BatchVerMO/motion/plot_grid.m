function plot_grid(xxsfyysf,overlap)
    xx_s=xxsfyysf{1};
    xx_f=xxsfyysf{2};
    yy_s=xxsfyysf{3};
    yy_f=xxsfyysf{4};
    xx_s_e=[]; xx_f_e=[]; yy_s_e=[]; yy_f_e=[];
    for i=1:length(xx_s)
        for j=1:length(yy_s)
            xx_s_e=[xx_s_e max(xx_s(i)-overlap(1),xx_s(1))];
            xx_f_e=[xx_f_e min(xx_f(i)+overlap(1),xx_f(end))];
            yy_s_e=[yy_s_e max(yy_s(j)-overlap(2),yy_s(1))];
            yy_f_e=[yy_f_e min(yy_f(j)+overlap(2),yy_f(end))];
        end
    end
    
    plot([0 d2], [xx_s_e'; xx_f_e']*[1 1], ':', 'color', .7*[1 1 1]);
    plot([yy_s_e'; yy_f_e']*[1 1],[0 d1], ':', 'color', .7*[1 1 1]);
