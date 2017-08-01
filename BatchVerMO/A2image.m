function image=A2image(A)
MfMAT=[];
A=A(:,any(A,1));
k = size(A,2);  C=ones(k,1);
Brainbow = A*C; Brainbow = reshape(Brainbow,300,400);
MfMAT=cat(3,MfMAT,Brainbow);
image=MfMAT;