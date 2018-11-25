init_fig;

% c
axes('position', [.05, .57, .95, .37]);
hold on;
plot(y, 'color', col{8}/255);
alpha(.7);
plot(true_c, 'color', col{3}/255, 'linewidth', 1.5);
plot(c_oasis, '-.', 'color', col{1}/255);
axis tight;
xlim([0, 2000]);
set(gca, 'xtick', [0, 25, 50, 75]*30);
set(gca, 'xticklabel', []);
set(gca, 'ytick', 0:2);
ylabel('Fluor.');
box off;
legend('Data', 'Truth', 'OASIS', 'location', 'northeast', 'orientation', 'horizental');
% s
axes('position', [.05, .18, .95, .37]);
hold on;
plot(true_s, 'color', col{3}/255, 'linewidth', 1.5);
plot(s_oasis, '-.', 'color', col{1}/255);
axis tight;
xlim([0, 2000]);
set(gca, 'xtick', [0, 25, 50, 75]*30);
set(gca, 'xticklabel', get(gca, 'xtick')/30);
set(gca, 'ytick', [0,1]);
xlabel('Time [s]');
ylabel('Activity.');