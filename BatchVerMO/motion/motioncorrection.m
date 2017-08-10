function [M_final,xxsfyysf]=motioncorrection(AsfromDaysCell,AsfromDaysPic,options_nonrigid,gridstartend)
%%
Y=AsfromDaysPic; 
clear AsfromDaysPic
AnumInFolder=numel(AsfromDaysCell);
sizes=cellfun(@(x) size(x,2),AsfromDaysCell);

M=cell(1,AnumInFolder);  M_central=cell(1,AnumInFolder); 
ind_del=cell(1,AnumInFolder);

% the only left out one;
M{1}{1}=AsfromDaysCell{1};  ind_del{1}{1}=false(1,sizes(1));
M_central{1}{1}=centralA(AsfromDaysCell{1});

for ia=2:AnumInFolder
    shifts_up_1=cell(1,ia-1); shifts_up_2=cell(1,ia-1);

    siz_oneday=sizes(ia);
    % The day's own data
    M{ia}{ia}=AsfromDaysCell{ia};
    M_central{ia}{ia}=centralA(AsfromDaysCell{ia});
    ind_del{ia}{ia}=false(1,siz_oneday);
    % Previous days
    M_buffer=[];
    for io=ia-1:-1:1
        if io==ia-1
            Y_template=Y(:,:,io);
        else
            Y_template=A2image(cat(2,M{io}{1:ia-1}),size(Y,1),size(Y,2));
        end
        if io==ia-1
            Y_toregister=Y(:,:,ia);
        else
            Y_toregister=A2image(cat(2,M{io+1}{io+1},M{io+1}{ia}),size(Y,1),size(Y,2));
        end
 
        %tic; [M{ia},shifts{ia},~,xxsfyysf,ind_del{ia}] = normcorre_BatchVer(Y_ex_oneday,options_nonrigid,Y_oneday,siz_ex_oneday,As_ex_oneday,startendgrid,update_num); toc
        tic; [M_temp,shifts,shifts_up_1{io},shifts_up_2{io},~,xxsfyysf,ind_del_temp] = ...
            normcorre_BatchVer(Y_toregister,options_nonrigid,Y_template,siz_oneday,M{io+1}{ia},gridstartend); toc
        M{io}{ia}=M_temp{1};
        M_central{io}{ia}=centralA(M_temp{1});
        ind_del{io}{ia}=ind_del_temp{1};

        % inverse consec
        shifts_up_1_tmp=sum(cat(3,shifts_up_1{io:ia-1}),3); shifts_up_2_tmp=sum(cat(3,shifts_up_2{io:ia-1}),3);
        Mf_temp=[];
        ind_del{ia}{io}=false(1,sizes(io));
        for ni=1:sizes(io)
            Y_one_neuron=reshape(AsfromDaysCell{io}(:,ni),size(Y,1),size(Y,2));
            Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)) = ...
                imwarp(Y_one_neuron(gridstartend(1):gridstartend(2),gridstartend(3):gridstartend(4),gridstartend(5):gridstartend(6)),-cat(3,shifts_up_2_tmp.*(-1),shifts_up_1_tmp.*(-1)),options_nonrigid.shifts_method);
            if any(Y_one_neuron)==0
                ind_del{ia}{io}(ni)=true;
            end
            Mf_temp=[Mf_temp reshape(Y_one_neuron,[],1)];
        end
        M{ia}{io}=Mf_temp;
        M_central{ia}{io}=centralA(M{ia}{io});
    end
end

%%
ind_del_full_cell=cellfun(@(x) cell2mat(x), ind_del, 'UniformOutput',0);
ind_del_full=sum(reshape(cell2mat(ind_del_full_cell),1,sum(sizes),[]),3)>0;    %sum(cat(3,ind_del{:}))>0;
ind_del_full_cell=mat2cell(~ind_del_full,1,sizes);
N_eachday=cellfun(@(x) sum(x), ind_del_full_cell);
M_concat=cellfun(@(x) cat(2,x{:}), M, 'UniformOutput',0);
M_del = cellfun(@(x) x(:,~ind_del_full), M_concat, 'UniformOutput',0);
M_final = cellfun(@(x) mat2cell(x, [size(x,1)], N_eachday), M_del, 'UniformOutput',0);