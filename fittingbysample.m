function [ff,datafitted]= fittingbysample (data,MZ,RT, ll)




ind=1;
outliers=[];
vector=[1:size(data,1)];
iter=0;


Cc=cell(size(data,1),1);
for k=1:size(data,1)
    Cc{k,1}=['error' int2str(k)];
end

Cc=cellstr(Cc);


while ~isempty(ind)
    iter=iter+1;
    datafitted=[];
    for i=1:size(data,2)



        [f,gof] = fit( [MZ(:,i),RT(:,i)], data(:,i), 'lowess' ,'Robust','LAR');
        ff{i,1}=f;
        datafitted(:,i)=f(MZ(:,i),RT(:,i));
    end
    
    SE= data-datafitted;

        [~, ~ ,~ ,indice] = testboxplot(SE, Cc, ll, 0,0)
        ind=indice(:,1);
    
%     [ind]=find(abs(SE(:))>3)
    [~,row]=data2vector(SE);
   
    
    outliers=[outliers vector(unique(row(ind)))];
    vector(unique(row(ind)))=[];
    
    data(unique(row(ind)),:)=[];
    MZ(unique(row(ind)),:)=[];
    RT(unique(row(ind)),:)=[];
end
% figure()
% hist(SE)

