function image=A2image(A,d1,d2)

A=A(:,any(A,1));
k = size(A,2);
sA = sum(A,1);
A = bsxfun(@rdivide, A, sA)*prctile(sA, 5);
C=ones(k,1);
Brainbow = 1-A*C; Brainbow = reshape(Brainbow,d1,d2);
image=Brainbow;