function c = fig2_demo_deconvolveAR1(y, g, lam, video, trueSpikes)
%% extract neural activity from a fluorescence trace using an active set
% method for sparse nonnegative deconvolution

%% inputs:
%   y:  T*1 vector, raw fluorescence trace
%   g:  scalar, AR(1) coefficient
%   lam:  scalar, tuning parameter for sparsity

%% outputs
%   c: T*1 vector, denoised trace

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% proted from the Python implementation from Johannes Friedrich

%% initialization
y = reshape(y, [], 1);
if ~exist('g', 'var') || isempty(g);   g = 0.95; end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end

len_active_set = length(y);
Aset = ones(len_active_set, 4); % each row is one pool: (vi, wi, t, l)
Aset(:,1) = y - lam;            % vi
Aset(:,3) = (1:len_active_set);  % ti
Aset(end, 1) = y(end) - lam/(1-g);
Aset(end, 3) = len_active_set;
Aset(1,1) = max(Aset(1,1), 0);

%% run OASIS
ii = 1;
counter = 1;
if video
    figure('papersize', [15, 6]);
    init_fig;
    avi_file = VideoWriter('fig/Video.avi');
    avi_file.FrameRate = 10;
    avi_file.open();
else
    figure('papersize', [3,3]);
    init_fig;
end
while ii < len_active_set
    % find the active set
    while (ii<len_active_set) && (Aset(ii+1,1)>=Aset(ii,1)*g^(Aset(ii,4)))
        ii = ii + 1;
    end
    
    if ii == len_active_set; break; end
    
    if video
        cb_video(y, Aset, g, ii, trueSpikes);
        avi_file.writeVideo(getframe(gcf));
    else
        cb(y, Aset, g, counter, ii, trueSpikes);
    end
    counter = counter+1;
    %% merge pools
    while ii>0 && (Aset(ii+1,1) < Aset(ii,1)*g^(Aset(ii,4)))
        temp = Aset(ii,2) + Aset(ii+1,2)*(g^(2*Aset(ii, 4)));
        Aset(ii,1) = (Aset(ii,1)*Aset(ii,2) + Aset(ii+1,1)*Aset(ii+1, 2)* ...
            (g^(Aset(ii,4)))) / temp;
        if ii==1
            Aset(ii,1) = max(Aset(ii,1), 0);
        end
        Aset(ii, 2) = temp;
        Aset(ii, 4) = Aset(ii, 4) + Aset(ii+1, 4);
        Aset(ii+1, :) = [];
        ii = ii - 1;
    end
    
    ii = ii + 1;
    len_active_set = size(Aset,1);
end

%% construct solution for all t
c = zeros(size(y));
c(Aset(:, 3)) = Aset(:, 1);
for ii=1:len_active_set
    t0 = Aset(ii,3);
    tau = Aset(ii, 4);
    c(t0:(t0+tau-1)) = Aset(ii,1) * (g.^(0:(tau-1)));
end
if video
    avi_file.close();
end
end

function cb(y, Aset, g, counter, current, trueSpikes)
%% function for taking snapshot of OASIS procedures
if ~exist('fig', 'dir')
    mkdir fig;
end
len_active_set = size(Aset,1);
% construct solution
c = zeros(size(y));
Aset(:,1) = max(0, Aset(:, 1));
c(Aset(:, 3)) = Aset(:, 1);
for ii=1:len_active_set
    t0 = Aset(ii,3);
    tau = Aset(ii, 4);
    c(t0:(t0+tau-1)) = Aset(ii,1) * (g.^(0:(tau-1)));
end
% plot results
clf; hold on;
% use different color to separate pools
ymax = 1.2;
for m=1:2:len_active_set
    t0 = Aset(m, 3);
    t1 = Aset(m, 3) + Aset(m, 4);
    fill([t0+0.1, t1, t1, t0+0.1], [0.01, 0.01, ymax, ymax], [1,1,1]*0.9, 'edgecolor', 'w');
end
temp = find(trueSpikes);
for m=1:length(temp);
    plot([1, 1]*temp(m), [0, ymax], 'r', 'linewidth', 3);
end
ind = 48-ceil((y-min(y))/range(y)*16);
col_map = jet();
col = col_map(ind, :);
% current estimation of c
plot(c, 'color', 'k', 'linewidth', 2);
scatter(1:length(y), c, 100, col, 'filled', 'markeredgecolor', 'k');
% scatter plot the current pool
t = Aset(current, 3)+(1:Aset(current, 4));
plot(t, c(t), '+b', 'markersize', 10);
ylim([-0, ymax]);
xlim([1, length(y)-1]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
box off;
pause(.1);
saveas(gcf, sprintf('fig/%d.pdf', counter));
end


function cb_video(y, Aset, g, current, trueSpikes)
%% function for taking snapshot of OASIS procedures and save the results as a video
if ~exist('fig', 'dir')
    mkdir fig;
end
len_active_set = size(Aset,1);
% construct solution
c = zeros(size(y));
c(Aset(:, 3)) = Aset(:, 1);
for ii=1:len_active_set
    t0 = Aset(ii,3);
    tau = Aset(ii, 4);
    c(t0:(t0+tau-1)) = Aset(ii,1) * (g.^(0:(tau-1)));
end
s = c;
for t=2:length(y)
    s(t) = c(t) - g*c(t-1);
end
% plot results
clf;
axes('position', [0.05, 0.55, 0.9, 0.4]); cla; hold on;
ymax = 1.8;
% use different color to separate pools
for m=1:2:len_active_set
    t0 = Aset(m, 3);
    t1 = Aset(m, 3) + Aset(m, 4);
    fill([t0, t1, t1, t0], [0.01, 0.01, ymax, ymax], [1,1,1]*0.9, 'edgecolor', 'w');
end
temp = find(trueSpikes);
for m=1:length(temp);
    plot([1, 1]*temp(m), [0, ymax], 'r', 'linewidth', 3);
end
ax = gca;
ax.XAxisLocation = 'origin';
plot([0, 0], [-0.6, 0], 'w', 'linewidth', 4);
plot([0, length(y)], [0, 0], 'k', 'linewidth', 4);
plot([0, 0], [0, ymax], 'k', 'linewidth', 4);

ind = 48-ceil((y-min(y))/range(y)*16);
col_map = jet();
col = col_map(ind, :);
% current estimation of c
plot(c, 'color', 'k', 'linewidth', 2);
scatter(1:length(y), c, 80, col, 'filled', 'markeredgecolor', 'k');
% scatter plot the current pool
t = Aset(current, 3)+(1:Aset(current, 4));
plot(t, c(t), '+b', 'markersize', 10);
ylim([-0.6, ymax]);
xlim([0, length(y)-1]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
ylabel('Fluorescence');

%% plot spike trains
axes('position', [0.05, 0.1, 0.9, 0.4]); cla; hold on;
ymax = 1.5;
% use different color to separate pools
for m=1:2:len_active_set
    t0 = Aset(m, 3);
    t1 = Aset(m, 3) + Aset(m, 4);
    fill([t0, t1, t1, t0], [0.01, 0.01, ymax, ymax], [1,1,1]*0.9, 'edgecolor', 'w');
end
temp = find(trueSpikes);
for m=1:length(temp);
    plot([1, 1]*temp(m), [0, ymax], 'r', 'linewidth', 3);
end
ax = gca;
ax.XAxisLocation = 'origin';
plot([0, 0], [-0.6, 0], 'w', 'linewidth', 4);
plot([0, length(y)], [0, 0], 'k', 'linewidth', 4);
plot([0, 0], [0, ymax], 'k', 'linewidth', 4);
s = c;
plot([1,1], [0, s(1)], 'k', 'linewidth', 3);
for t=2:length(y)
    s(t) = c(t) - g*c(t-1);
    plot([t, t], [0, s(t)], 'k', 'linewidth', 3);
end
scatter(1:length(y), s, 80, col, 'filled', 'markeredgecolor', 'k');
% scatter plot the current pool
t = Aset(current, 3)+(1:Aset(current, 4));
plot(t, s(t), '+b', 'markersize', 10);
ylim([-0.6, ymax]);
xlim([0, length(y)-1]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
ylabel('Spikes');
xlabel('Time');
end