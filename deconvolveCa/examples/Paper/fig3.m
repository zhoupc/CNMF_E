%%
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]};

%% AR(1)
g = .95;
[Y, trueC, trueSpikes] = gen_data();
N = size(Y, 1); 

% run deconvolution
lam = 2.4;
[c_oasis, s_oasis] = oasisAR1(Y(1,:), g, lam);
[c_foopsi, s_foopsi] = foopsi(Y(1,:), g, lam);

% plot results
figure('papersize', [15, 4]);
init_fig;

% c
axes('position', [.05, .57, .95, .37]);
hold on;
plot(Y(1,:), 'color', col{8}/255);
alpha(.7);
plot(trueC(1,:), 'color', col{3}/255, 'linewidth', 1.5);
plot(c_oasis, 'color', col{1}/255);
plot(c_foopsi, '--', 'color', col{7}/255);
axis tight;
xlim([0, 2000]);
set(gca, 'xtick', [0, 25, 50, 75]*30);
set(gca, 'xticklabel', []);
set(gca, 'ytick', 0:2);
ylabel('Fluor.');
box off;
legend('Data', 'Truth', 'OASIS', 'CVX', 'location', 'northeast', 'orientation', 'horizental');
% s
axes('position', [.05, .18, .95, .37]);
hold on;
plot(trueSpikes(1,:), 'color', col{3}/255, 'linewidth', 1.5);
plot(s_oasis, 'color', col{1}/255);
plot(s_foopsi, '--', 'color', col{7}/255);
axis tight;
xlim([0, 2000]);
set(gca, 'xtick', [0, 25, 50, 75]*30);
set(gca, 'xticklabel', get(gca, 'xtick')/30);
set(gca, 'ytick', [0,1]);
xlabel('Time [s]');
ylabel('Activity.');

saveas(gcf, 'fig/traceAR1.pdf');

%% time it
fprintf('\n**************AR 1**************\n'); 
tic;
for m=1:N
    [c_oasis, s_oasis] = oasisAR1(Y(m,:), g, lam);
end
fprintf('OASIS: %.3f seconds\n', toc);

tic;
for m=1:N
    [c_foopsi, s_foopsi] = foopsi(Y(1,:), g, lam); %#ok<ASGLU>
end
fprintf('FOOPSI: %.3f seconds\n', toc);
fprintf('\n**************AR 1**************\n'); 

%% AR(2)
% simulation
g = [1.7, -0.712];
sn = 1;
seed = 3;
[Y, trueC, trueSpikes] = gen_data(g, sn, [], [], [], [], [], seed);
[N, T] = size(Y);

% deconvolution
lam = 25;
[c_onnls, s_onnls] = onnls(Y(1,:), g, lam);
[c_foopsi, s_foopsi] = foopsi(Y(1,:), g, lam);
figure('papersize', [15, 4]);
init_fig;

% plot results
figure('papersize', [15, 4]);
init_fig;

% c
axes('position', [.05, .57, .46, .37]);
hold on;
plot(Y(1,:), 'color', col{8}/255, 'linewidth', 0.5);
alpha(.7);
plot(trueC(1,:), 'color', col{3}/255, 'linewidth', 1.5);
plot(c_onnls, 'color', col{1}/255);
plot(c_foopsi, '--', 'color', col{7}/255);
axis tight;
xlim([0, 1200]);
set(gca, 'xtick', 0:300:1500);
set(gca, 'xticklabel', []);
ylim(round([1+min(Y(1,:)), max(Y(1,:))-0.5]));
ylabel('Fluor.');
box off;
% s
axes('position', [.05, .18, .46, .37]);
hold on;
plot(trueSpikes(1,:), 'color', col{3}/255, 'linewidth', 1.5);
plot(s_onnls, 'color', col{1}/255);
plot(s_foopsi, '--', 'color', col{7}/255);
axis tight;
xlim([0, 1200]);
set(gca, 'xtick', 0:300:1500);
set(gca, 'xticklabel', get(gca, 'xtick')/30);
set(gca, 'ytick', [0,1]);
xlabel('Time [s]');
ylabel('Activity.');

%% timeit
fprintf('\n**************AR 2**************\n'); 
tic;
for m=1:N
    [c_onnls, s_onnls] = onnls(Y(m,:), g, lam);
end
fprintf('online NNLS: %.3f seconds\n', toc);

tic;
for m=1:N
    [c_foopsi, s_foopsi] = foopsi(Y(m,:), g, lam);
end
fprintf('FOOPSI: %.3f seconds\n', toc);
fprintf('\n**************AR 2**************\n'); 

















