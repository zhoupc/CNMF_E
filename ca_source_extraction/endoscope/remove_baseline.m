function [y, b] = remove_baseline(y)
% estiamte baseline of the calcium traces and subtract it 
sz = size(y); 
y = reshape(y, 1, []); 
noise = get_noise_fft(reshape(y, 1, [])); 
y_diff = [-1, diff(y)]; 
b = median(y(and(y_diff>=0, y_diff<noise))); 
y = reshape(y-b, sz); 