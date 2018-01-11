%% parameters for simulation 
gam = [1.58, -0.6]; 
noise = .5; 
T = 455; 
framerate = 30; 
firerate = 2; 
seed = 3; 

%% simulate calcium fluorescience 
[Y, trueC, trueSpikes] = gen_data(gam, noise, T, framerate, ...
    firerate, [], [], seed); 

%% plot results 
figure('papersize', [15, 2.5]); 
init_fig; 
hold on; 
col = {[0, 114, 176]; ...
    [0, 158, 115]; ...
    [213, 94, 0]}; 

% fluorescence trace 
plot(1:T, Y(1,:)/3, 'o', 'color',  uint8(col{2}));

% calcium trace 
plot(1:T, trueC(1,:)/3, 'color', 'k'); % uint8(col{1})); 
% spike train 
tsp = find(trueSpikes(1, :)); 
for m=1:length(tsp)
    plot([1,1]*tsp(m), [0, 1], 'color', uint8(col{3})); 
end
axis tight; 
xlabel('Time'); 
ylabel('Fluorescence'); 
legend('y', 'c', 's'); 
set(gco, 'fontweigth', 'bold'); 
saveas(gcf, 'fig/model.pdf'); 
