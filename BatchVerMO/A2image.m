function image=A2image(A,d1,d2,textOrNot,color,PoporNOt,handle)

if nargin<4
    textOrNot=false;
end

if nargin<5
    color=[];
end

if nargin<6
    PoporNOt=false;
end

A=A(:,any(A,1));
k = size(A,2);
sA = sum(A,1);
A = bsxfun(@rdivide, A, sA)*prctile(sA, 5);
C=ones(k,1);

if isempty(color)==false
    if or(strcmp(color,'magenta'),strcmp(color,'red'))
        color_palet = 1-[1 0 1];
    elseif strcmp(color,'green')
        color_palet = 1-[0 1 0];
    end
        nColors = repmat(color_palet,size(A,2),1); 
        Ir = diag(nColors(:,1));
        Ig = diag(nColors(:,2));
        Ib = diag(nColors(:,3));
    
    Brainbow = 1-cat(3, A*Ir*C, A*Ig*C, A*Ib*C);
    Brainbow = reshape(Brainbow,d1,d2,3);
else
    Brainbow = 1-A*C;
    Brainbow(Brainbow<0)=0;
    Brainbow = reshape(Brainbow,d1,d2);
    %image=imshow(Brainbow); 
end
image=Brainbow;
if PoporNOt==true
    if nargin<7
        figure
    else
        axes(handle)
    end
    imshow(Brainbow)
end

if textOrNot==true
    Atemp=reshape(A,d1,d2,k);
    Position=zeros(2,k);
    for i=1:k
        Atemp_=Atemp(:,:,i);
        [row_ind,col_ind] = find(Atemp_>0);
        Position(2,i)=mean(row_ind);
        Position(1,i)=mean(col_ind);
    end
    text(Position(1,:),Position(2,:),cellstr(num2str((1:k)'))','Color','black')
    F = getframe(gca);
    image=F.cdata; 
end
