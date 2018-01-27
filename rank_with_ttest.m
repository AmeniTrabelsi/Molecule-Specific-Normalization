function [ranks]= rank_with_ttest(data1,data2)

[~,p]=ttest2(data1',data2');
[~,ind]=sort(p,'ascend');
[~,pos]=sort(ind,'ascend');

ranks(:,1)=1:size(data1,1);
ranks(:,2)=pos;
