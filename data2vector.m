function [v,r,c]=data2vector(data)

dim1=size(data,1);
dim2=size(data,2);

r=repmat([1:dim1]',dim2,1);
c=repmat([1:dim2],dim1,1);
c=c(:);
v=data(:);