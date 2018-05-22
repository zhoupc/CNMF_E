function Position=CalText_AddText(Position,A,d1,d2)
    if nargin<2
        k=size(Position,2);
        text(Position(1,:),Position(2,:),cellstr(num2str((1:k)'))','Color','black')
    else
        plotOrNot=Position;
        k=size(A,2);
        Atemp=reshape(A,d1,d2,k);
        Position=zeros(2,k);
        for i=1:k
            Atemp_=Atemp(:,:,i);
            [row_ind,col_ind] = find(Atemp_>0);
            Position(2,i)=mean(row_ind);
            Position(1,i)=mean(col_ind);
        end
        if plotOrNot
            text(Position(1,:),Position(2,:),cellstr(num2str((1:k)'))','Color','black');
        end
    end