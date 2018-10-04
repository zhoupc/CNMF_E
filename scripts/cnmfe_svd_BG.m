clear Ysignal;
tic; 
Ybg = Y-neuron.A*neuron.C;
b0 = mean(Ybg, 2); 
[u, s, v] = svdsecon(bsxfun(@minus, Ybg, b0), nb); 
Ybg = bsxfun(@plus, u*s*v', b0); 
Ysignal = Y - Ybg;
