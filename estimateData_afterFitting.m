function [data]=estimateData_afterFitting(ff,MZ,RT)

for i=1:size(ff,1)
    f=ff{i,1};
    data(:,i)=f(MZ(:,i),RT(:,i));
end