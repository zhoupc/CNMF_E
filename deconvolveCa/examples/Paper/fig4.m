%%
close all; clc; 
addpath('./scripts'); 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]};
figure('papersize', [14.4, 7.2]); 
init_fig; 

%% load the test data  
gam = 0.95; 
sn = .3;
T = 1500; 
lam = 0; 
load fig4_data; 
y = reshape(y, [],1); 
g = .9; 
g0 = g; 
ax0 = .06; 
ax1 = .2; 

%% plot initial result 
[solution, spks, active_set] = oasisAR1(y, g, lam); 
axes('position', [ax1, 0.87, 1-ax1, .12]); 
fig4_plot_trace; 
legend('Data', 'Truth', 'Estimate', 'orientation', 'horizental',...
    'location', [0.01+ax1, 0.96, 0.35, 0.0345]); 

%% solve for lambda 
c = solution; 
temp = zeros(size(c)); 
len_active_set = size(active_set,1); 
for ii=1:len_active_set
    ti = active_set(ii, 3); 
    li = active_set(ii, 4); 
    idx = 0:(li-1); 
    temp(ti+idx) = (1-g^li)/ active_set(ii,2) * g.^(idx);   
end
res = y - solution; 
aa = temp'*temp; 
bb = res'*temp; 
cc = res'*res - sn^2*T; 
active_set_0 = active_set; 
ll = (-bb + sqrt(bb^2-aa*cc)) / aa; 
active_set_0(:,1) = active_set_0(:,1) - ll*(1-g.^(active_set_0(:,4)));
spks = 0*spks; 
for ii=1:len_active_set
    ti = active_set(ii, 3); 
    li = active_set(ii, 4); 
    idx = 0:(li-1); 
    solution(ti+idx) = active_set_0(ii,1)/active_set_0(ii,2) * (g.^(idx));  
    if ii>1
        spks(ti) = solution(ti) - g*solution(ti-1); 
    end
end

% plot results after updating lambda, but before running oasis to fix
% violations 
axes('position', [ax1, .73, 1-ax1, .12]); 
fig4_plot_trace; 

% plot the depence of RSS on lambda 
lam_vec = linspace(0, 2.0, 100); 
RSS_vec = zeros(size(lam_vec)); 
for m=1:length(lam_vec); 
    RSS_vec(m) = norm((res-lam_vec(m)*temp),2)^2; 
end

axes('position', [ax0, .75, .08, .12]); hold on; 
plot(lam_vec, RSS_vec, 'color', col{2}/255); 

% the optimal value 
thresh = sn.^2 * T; 
plot([-.1, ll], thresh*[1,1], 'k'); 
plot([ll, ll], [100, thresh], 'k'); 
scatter(0, RSS_vec(1), 50, col{1}/255); 
xlim([0, 1.7]); 
ylim([114, 138]); 
set(gca,'ytick', thresh); 
set(gca, 'yticklabel', '\sigma^2T'); 
set(gca, 'xtick', [0,lam+ll]); 
set(gca, 'xticklabel', {0, '\lambda^*'}); 

%% plot result after rerunning oasis to fix violations
lam = lam + ll; 
[solution, spks, active_set] = oasisAR1(y, g, lam, [], active_set_0); 
axes('position', [ax1, .59, 1-ax1, .12]); 
fig4_plot_trace; 
ylabel('Fluorescence'); 

%% solve for gamma 
[~, ~, g] = update_g(y, active_set,g, lam); 
axes('position', [ax0, .47, .08, .12]); hold on; 
g_vec = linspace(.85, 0.99); 
rss_vec = compute_rss_g(g_vec, y, active_set, lam); 
plot(g_vec, rss_vec, 'color', col{2}/255); 
xlim([min(g_vec), max(g_vec)]); 
scatter(g0, compute_rss_g(g0, y, active_set, lam), 50, col{1}/255, 'filled'); 
set(gca, 'xtick', g); 
set(gca, 'xticklabel', '\gamma^*'); 
set(gca,'ytick', thresh); 
set(gca, 'yticklabel', '\sigma^2T'); 
% plot result after updating gamma, but before rerunning oasis to fix
% violations 
len_active_set = size(active_set, 1);
ma = max(active_set(:, 4)); 
h = g.^((0:ma)'); 
hh = cumsum(h.*h); 
for ii=1:len_active_set
    ti = active_set(ii,3); 
    li = active_set(ii,4); 
    idx = ti:(ti+li-1); 
    active_set(ii,1) = (y(idx)-lam*(1-g))'*h(1:li); 
    active_set(ii,2) = hh(li); 
end
for ii=1:length(active_set_0)
    ti = active_set_0(ii,3); 
    li = active_set_0(ii,4); 
    idx = ti:(ti+li-1); 
    solution(idx) = active_set_0(ii,1)/active_set_0(ii,2) * h(1:li); 
end
solution(solution<0) = 0; 
spks = [0; solution(2:end)-g*solution(1:(end-1))]; 
axes('position', [ax1, .45, 1 - ax1, .12]); hold on; 
fig4_plot_trace; 

%% plot result after rerunning oasis to fix violatiosn 
[solution, spks, active_set] = oasisAR1(y, g, lam,[], active_set); 
axes('position', [ax1, .31, 1 - ax1, .12]); hold on; 
fig4_plot_trace; 

%% do more iterations 
for miter=1:3
    [~, active_set, lam, ~] = update_lam(y, solution, active_set, g, lam, thresh); 
    [solution, active_set, g, spks] = update_g(y, active_set, g, lam); 
end
axes('position', [ax1, .07,  1 - ax1, .12]); 
hold on; 
fig4_plot_trace; 
sol_given_g = constrained_oasisAR1(y, .95, sn); 
estimated_g = estimate_parameters(y, 1); 
fprintf('\n*******************\n'); 
fprintf('estimated gamma via autocorrelation: %.4f\n', estimated_g); 
fprintf('optimized gamma                    : %.4f\n', g); 
fprintf('\n*******************\n'); 
sol_PSD_g = constrained_oasisAR1(y, estimated_g, 0); 
h1 = plot(sol_given_g, '--', 'color', col{7}/255); 
h2 = plot(sol_PSD_g, 'color', col{6}/255); 
legend([h1, h2], 'true \gamma', '\gamma from autocovariance', ...
    'orientation', 'horizontal', 'location', [0.01+ax1, 0.16, 0.35, 0.0345]); 
saveas(gcf, 'fig/opt_g+l_new.pdf')