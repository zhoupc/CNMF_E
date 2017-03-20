%%
clear; clc; close all; 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]};

%% AR (1) 
% simulation 
g = .95; 
sn = .3; 
[Y, trueC, trueSpikes] = gen_data(g, sn, [], [],[], [], 1, 10); 
y = Y(1, :); 

% run deconvolution 
[c, s] = constrained_foopsi_cvx(y, g, sn); 
[c_t, s_t] = oasisAR1(y, g, 0, .55); 

% check the dependence on smin 
smin = 0:0.1:1.1; 
res = zeros(length(smin), length(y)); 
for m=1:length(smin)
    [~, res(m, :)] = oasisAR1(y, g, 0, smin(m)); 
end
% plot results 
figure('papersize', [10,9]); 
init_fig; 
ymax = ceil(max(y)*2)/2; 
ymin = quantile(y(1:500), 0.02); 
% c
axes('position', [.13, .7, .86, .29]); hold on;  
plot(y, 'color', col{8}/255); 
alpha(.7); 
plot(trueC(1, :), 'color', col{3}/255, 'linewidth', 3); 
plot(c, 'color', col{1}/255); 
plot(c_t, 'color', col{2}/255); 
legend('Data', 'Truth', 'Thresh', 'L1', 'orientation', 'horizental'); 
ylim([ymin, ymax]); 
xlim([0, 452]); 
set(gca, 'xtick', 0:150:500); 
set(gca, 'xticklabel', []); 
set(gca, 'ytick', 0:1:ymax); 
ylabel('Fluor.'); 

%s 
axes('position', [.13, .39, .86, .29]); hold on;  
for m=find(trueSpikes(1,1:500))
    plot([m, m], [0, 1], 'color', col{3}/255); 
end
plot([0, 450], [0, 0], 'color', col{3}/255); 
for m=find(s(1:500))
    plot([m, m], [0, s(m)]+2.5, 'color', col{1}/255); 
end
plot([0, 450], [1, 1]*2.5, 'color', col{1}/255); 
for m=find(s_t(1:500)')
    plot([m, m], [0, s_t(m)]+1.25, 'color', col{2}/255); 
end
plot([0, 450], [1, 1]*1.25, 'color', col{2}/255); 

ylim([0, 3.5]); 
xlim([0, 452]); 
set(gca, 'xtick', 0:150:500); 
set(gca, 'xticklabel', []); 
set(gca, 'ytick', 0:1.25: 2.5);
set(gca, 'yticklabel', {'Thresh.', 'Truth', 'L1'})

% dependence on smin 
axes('position', [.13, .08, .86, .29]); hold on;  
for m=find(trueSpikes(1,1:500))
    plot([m, m], [-.08, -0.16], 'r'); 
end
for rr=1:length(smin)
    for m=find(res(rr, 1:500))
        plot([m, m], rr*.1+[-0.14, -0.06], 'color', 'k');
    end
end

ylim([-0.2, smin(end)]); 
xlim([0, 452]); 
set(gca, 'xtick', 0:150:500); 
set(gca, 'xticklabel', get(gca, 'xtick')/30)
set(gca, 'ytick', 0:0.5:1); 
ylabel('s_{min}'); 

saveas(gcf, 'fig/threshAR1.pdf'); 


%% AR (2) 
% simulation 
g = [1.7, -0.712]; 
sn = 1; 
seed = 1; 
[Y, trueC, trueSpikes] = gen_data(g, sn, [],[],[],[], [], seed); 
% Y(:, 1:150) = []; 
% trueC(:,1:150) = []; 
% trueSpikes(:,1:150) = []; 
y = Y(1, :); 

% run deconvolution 
[c, s] = constrained_foopsi_cvx(y, g, sn); 
[c_t, s_t] = oasisAR2(y, g, 0, .55); 

% check the dependence on smin 
smin = 0:0.1:1.1; 
res = zeros(length(smin), length(y)); 
for m=1:length(smin)
    [~, res(m, :)] = oasisAR2(y, g, 0, smin(m)); 
end
% plot results 
figure('papersize', [10,9]); 
init_fig; 
ymax = ceil(max(y)*2)/2; 
ymin = quantile(y(1:500), 0.02); 
% c
axes('position', [.13, .7, .86, .29]); hold on;  
plot(y, 'color', col{8}/255); 
alpha(.7); 
plot(trueC(1, :), 'color', col{3}/255, 'linewidth', 3); 
plot(c, 'color', col{1}/255); 
plot(c_t, 'color', col{2}/255); 
legend('Data', 'Truth', 'Thresh', 'L1', 'orientation', 'horizental'); 
ylim([ymin, ymax]); 
xlim([0, 452]); 
set(gca, 'xtick', 0:150:500); 
set(gca, 'xticklabel', []); 
set(gca, 'ytick', 0:1:ymax); 
ylabel('Fluor.'); 

%s 
axes('position', [.13, .39, .86, .29]); hold on;  
for m=find(trueSpikes(1,1:500))
    plot([m, m], [0, 1], 'color', col{3}/255); 
end
plot([0, 450], [0, 0], 'color', col{3}/255); 
for m=find(s(1:500))
    plot([m, m], [0, s(m)]+2.5, 'color', col{1}/255); 
end
plot([0, 450], [1, 1]*2.5, 'color', col{1}/255); 
for m=find(s_t(1:500)')
    plot([m, m], [0, s_t(m)]+1.25, 'color', col{2}/255); 
end
plot([0, 450], [1, 1]*1.25, 'color', col{2}/255); 

ylim([0, 3.5]); 
xlim([0, 452]); 
set(gca, 'xtick', 0:150:500); 
set(gca, 'xticklabel', []); 
set(gca, 'ytick', 0:1.25: 2.5);
set(gca, 'yticklabel', {'Thresh.', 'Truth', 'L1'})

% dependence on smin 
axes('position', [.13, .08, .86, .29]); hold on;  
for m=find(trueSpikes(1,1:500))
    plot([m, m], [-.08, -0.16], 'r'); 
end
for rr=1:length(smin)
    for m=find(res(rr, 1:500))
        plot([m, m], rr*.1+[-0.14, -0.06], 'color', 'k');
    end
end

ylim([-0.2, smin(end)]); 
xlim([0, 452]); 
set(gca, 'xtick', 0:150:500); 
set(gca, 'xticklabel', get(gca, 'xtick')/30)
set(gca, 'ytick', 0:0.5:1); 
ylabel('s_{min}'); 

saveas(gcf, 'fig/threshAR2.pdf'); 
























