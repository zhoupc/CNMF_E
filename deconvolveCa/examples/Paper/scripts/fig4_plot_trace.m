plot(y, 'color', col{8}/255, 'linewidth', 1); 
hold on; 
plot(trueC, 'color', col{3}/255); 
plot(solution, 'color', col{1}/255); 
text(700, 2.8, sprintf('Correlation: %.3f', corr(trueS, spks)), ...
    'fontweight', 'bold', 'fontsize', 16); 
xlim([0, 1500]); 
ylim([min(y), 3.3]); 
ax = gca; 
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca, 'ytick', 0:2:max(y)); 
set(gca, 'xtick', 0:300:1200); 
set(gca, 'xticklabel', []); 
set(gca, 'yticklabel', [0,2]); 
box off; 