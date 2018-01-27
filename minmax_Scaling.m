function [normdata, minn, maxx]= minmax_Scaling(data)

normdata=data;

for i=1:size(data,2)
    SortedData = sort(data(:,i), 'ascend');
    TrimmedData = SortedData(floor(0.05*size(SortedData)):floor(size(SortedData)));
%     TrimmedData=normdata(:,i);
    minn(i)=min(TrimmedData);
    maxx(i)=max(TrimmedData);
    
    normdata(:,i)= (normdata(:,i)-minn(i))/(maxx(i)-minn(i));
end