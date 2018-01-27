function [rank]=rank_with_plsda( Data1,Data2)
addpath('D:\Users\ameni\Masters work\SurfaceFittingNormalization\PLSDA')
% addpath('D:\Users\ameni\Masters work\SurfaceFittingNormalization\libPLS_1.95')
% LDA=plslda([Data1 Data2]',[ones(size(Data1,2),1);2*ones(size(Data2,2),1)]);
% % rank(:,1)=[1:size(Data,1)];
% trah=LDA.tpLoadings;
% trah=abs(trah(1:size(trah,1)-1));

[plspara]=getplsda([Data1 Data2],[zeros(1,10) 5*ones(1,10)])
trah=plspara.vip;
[a,b]=sort(trah(:,1),'descend');
b(:,2)=[1:size(b,1)];
rank=sortrows(b);