function [data]= vector2data(vector, row)

for i=1:max(row)
    data(i,:)=vector(find(row==i));
end