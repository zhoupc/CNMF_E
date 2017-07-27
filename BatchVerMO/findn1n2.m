function [ns_storage] = findn1n2(A1,A2,correlation_thresh,max2max2nd,skewnessthresh)
% [ns_storage] = findn1n2(A1,A2,correlation_thresh,max2max2nd,skewnessthresh)
% Some input explained: only those neurons that have coor coef with other neurons with skewness over
%       skewnessthresh will be considered. See matlab function skewness().
%       The idea is that if this neuron is very independent, then its coor coef should only have very few big values.
%       While those skewness with a lot of overlapping with other neurons
%       will have many median level coor coef. skewness can very
%       efficiently distinguish these.
% Output: if there is no tracked neurons over days. The output is [].
% This function is the base for TrackingNeuronOverDays.m and Over_Days_findAnn.m

% by Shijie Gu, techel@live.cn, ShanghaiTech University, Harvard-MIT.

if isempty(correlation_thresh)
    correlation_thresh=0.8;
else
    correlation_thresh=correlation_thresh;
end

if isempty(max2max2nd)
    max2max2nd=1.1;
else
    max2max2nd=max2max2nd;
end

if isempty(skewnessthresh)
    skewnessthresh=5;
else
    skewnessthresh=skewnessthresh;
end


%A1=neuron1.A;  
A1=A1;
N1=1:size(A1,2);
%A2=neuron2.A;
A2=A2;
N2=1:size(A2,2);

%%  computing coor
% normalize
    %CoorMatrix=(A1)'*(A2)
    % A1norm=arrayfun(@(idx) norm(A1(:,idx)), 1:size(A1,2));
    % A2norm=arrayfun(@(idx) norm(A2(:,idx)), 1:size(A2,2));
    % CoefMatrix=bsxfun(@rdivide,CoorMatrix,A1norm');
    % CoefMatrix=bsxfun(@rdivide,CoefMatrix,A2norm);

temp1 = bsxfun(@times, A1>0, 1./sqrt(sum(A1>0)));
temp2 = bsxfun(@times, A2>0, 1./sqrt(sum(A2>0)));
CoefMatrix = temp1'*temp2;

%% filtering Coef
% (1) Each row/column must only have one value, its max.
%     There might be situations where row/col has the same max value though
%     rare. But this situation will be filtered out in the next step.
Filter1r=bsxfun(@eq,CoefMatrix,max(CoefMatrix,[],2));
Filter1c=bsxfun(@eq,CoefMatrix,max(CoefMatrix,[],1));
Filter1=and(Filter1r,Filter1c);
    % The following is for the situation when there are two same max.
    % Finding out rows with two maximums.
% [r,c]=find(Filter1);
% rc=[r,c];
% [r_sorted,r_sort_index]=sort(r,'ascend');
% rc_sorted=rc(r_sort_index,:);
% rowind=diff([0;r_sorted])==0;
% row_two_max=r_sorted(rowind);
% if length(row_two_max)>=1
%     protected=[];
%     protected2=[];
%         % (2) Assign one of the two max.
%     for i=1:length(row_two_max)
%         colind=find(rowind,i);
%         colind=colind(end);
%         col_ind=rc_sorted(colind-1,2);
%         col_ind_2=rc_sorted(colind,2);
%         if ~ismember(col_ind,protected)
%             Filter1(row_two_max(i),col_ind_2)=0;
%             protected=[protected col_ind]
%         else
%             if ~ismember(col_ind_2,protected2)
%                Filter1(row_two_max(i),col_ind)=0;
%                protected2=[protected2 col_ind_2]
%             else
%                Filter1(row_two_max(i),col_ind)=0;
%                Filter1(row_two_max(i),col_ind_2)=0;
%             end
%         end
%     end
% end

% (2) Confidence that the neuron is the neuron we trace comes from the fact that
%     (i) this correlation is unique.
%     (ii) it is spatitally not overlapped with others.

Coef_sorted=sort(CoefMatrix,2,'descend');
try
    ConfidentRow=(Coef_sorted(:,1)./Coef_sorted(:,2))>max2max2nd;
    Filter3r=repmat(ConfidentRow,1,size(A2,2));
catch
    Filter3r=true(size(Coef_sorted,1),1);
end

Coef_sorted2=sort(CoefMatrix,1,'descend');
try
    ConfidentCol=(Coef_sorted2(1,:)./Coef_sorted2(2,:))>max2max2nd;
    Filter3c=repmat(ConfidentCol,size(A1,2),1);
catch
    Filter3c=true(1,size(Coef_sorted,2));
end
Filter3=Filter3r&Filter3c;

%[skewness_sorted,ind]=sort(skewness(CoefMatrix,1,2),'descend');
skew_ind=skewness(CoefMatrix,1,2);
ConfidentRow_skew=skew_ind>=skewnessthresh;
skew_ind_Col=skewness(CoefMatrix,1,1);
ConfidentCol_skew=skew_ind_Col>=skewnessthresh;
Filter4r=repmat(ConfidentRow_skew,1,size(A2,2));
Filter4c=repmat(ConfidentCol_skew,size(A1,2),1);
Filter4=Filter4r&Filter4c;

Filter2=CoefMatrix>correlation_thresh;


FilterCoef_temp=and(Filter1,Filter2);
FilterCoef_temp=and(FilterCoef_temp,Filter3);
FilterCoef=and(FilterCoef_temp,Filter4);

% Calculating p-value
% t=CoefMatrix./sqrt(1-CoefMatrix.^2)*sqrt(size(A1,1)-2);
%(http://support.minitab.com/en-us/minitab-express/1/help-and-how-to/modeling-statistics/regression/how-to/correlation/methods-and-formulas/)

    % Let n be your sample size (size(A1,1))
    % Let v be your degrees of freedom (2)
    % Then: (https://stackoverflow.com/questions/10617050/how-to-calculate-p-value-for-t-test-in-matlab)
% pvalues = 2*(1-tcdf(abs(t),size(A1,1)-2));
% Filterpvalue = pvalues<=0.5;
%% 
% Filter= and(Filterpvalue,FilterCoef);
Filter=FilterCoef;
corrnumber=sum(sum(Filter,1),2);
fprintf('There are %.0f pairs of neurons tracked over two days.\n', corrnumber);

[n1,n2] = find(Filter);
ns_storage=[n1,n2];
% n1_the_rest=setdiff(N1,n1);% returns the data in N1 that is not in n1
% n2_the_rest=setdiff(N2,n2);
% 
% neuron1.APermuted=A1(:,[n1' n1_the_rest]);
% neuron2.APermuted=A2(:,[n2' n2_the_rest]);
end