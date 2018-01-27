function [normdata,standev,means]=Pareto_Scaling(data,dim)

normdata=data;
    if dim==2
     for i=1:size(data,2)   
         SortedData = sort(data(:,i), 'ascend');
         TrimmedData = SortedData(1:size(SortedData));
         means(i)=mean(TrimmedData);
         standev(i)=std(TrimmedData);
         normdata(:,i)=(normdata(:,i)-means(i))/sqrt(standev(i));
%          Bins = [min(data(:,i)):(max(data(:,i))-min(data(:,i)))/50: max(data(:,i))]
%          figure(i); plot(Bins, hist(data(:,i), Bins))
     end
     
    else
        for i=1:size(data,1)   
         means=mean(data(i,:));
         standev=std(data(i,:));
         normdata(i,:)=(normdata(i,:)-means)/sqrt(standev);
        end
    end
end