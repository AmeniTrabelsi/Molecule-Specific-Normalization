function [finaldata,normdata]= Surf_Fit(Data1,Data2, data_MZ, data_RT, thr)

% Input:
% Data1 and Data2: data to normalize.
% data_MZ: M/Z data corresponding to the data.
% data_RT: Retention Time data corresponding to the data.
% thr: Boxplot Threshold (ex: 1.5, 2...)
% 
% Output:
% 
% finaldata: data after normalization using MSN.


Data=[Data1, Data2];
[feature_ranks]=rank_with_ttest(Data1,Data2);
% Biomarkers are the first 27 metabolites in our data
BiomarkersSel=feature_ranks(find(feature_ranks(:,2)<=27),1);

% or 
% load('feature_ranks.mat')
% BiomarkersSel=feature_ranks(find(feature_ranks(:,2)<=27),1);
Housekeeping=setxor(1:size(Data,1),BiomarkersSel);
data_MZ_HK=data_MZ(Housekeeping,:);
data_RT_HK=data_RT(Housekeeping,:);

%% Scaling Data 


%MinMax Scaling
% [normdata,minn,maxx]=minmax_Scaling(Data,0);
[normdata,standev,means]=Pareto_Scaling(Data,2);


%%  Dividing by Median
medians=[];
testmedian=[];
testmediann=[];
normdataa=normdata+max(median(normdata,2));
medians=median(normdataa,2);
medians=abs(medians);
test=repmat(medians,1,size(Data,2));
testmediann=(normdataa)./(test);
testmedian=testmediann(Housekeeping,:);


%% Surface Fitting
parfor i =1:size(testmedian,2)
    [ff]=fittingbysample(testmedian(:,i),data_MZ_HK(:,i),data_RT_HK(:,i),thr);
    fff{i,1}=ff{1,1};

end


[finaldataN]=estimateData_afterFitting(fff,data_MZ,data_RT);



%% get Final Data
finaldata=(normdataa)./abs(finaldataN);

