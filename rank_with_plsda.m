function [rank]=rank_with_plsda( Data1,Data2)

addpath('./libPLS_1.95')
LDA=plslda([Data1 Data2]',[ones(size(Data1,2),1);2*ones(size(Data2,2),1)]);
trah=LDA.VIP;
trah=trah';


% or
% addpath('./PLSDA')
% [plspara]=getplsda([Data1 Data2],[zeros(1,10) 5*ones(1,10)])
% trah=plspara.vip;

[a,b]=sort(trah(:,1),'descend');
b(:,2)=[1:size(b,1)];
rank=sortrows(b);