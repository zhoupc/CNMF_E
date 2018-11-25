addpath('./scripts'); 
%% figures 
gam = .8; 
T = 25; 
noise = .2; 
seed = 3; 
[y, ~, trueSpikes] = gen_data(gam, noise, T, 1, 0.1, [], 1, seed); 

fig2_demo_deconvolveAR1(y, gam, 0.4, false, trueSpikes); 

%% video 
gam = .8; 
T = 130; 
noise = .2; 
seed = 3; 
[y, truth, trueSpikes] = gen_data(gam, noise, T, 1, 0.1, [], 1, seed); 
fig2_demo_deconvolveAR1(y, gam, 0.4, true, trueSpikes); 
