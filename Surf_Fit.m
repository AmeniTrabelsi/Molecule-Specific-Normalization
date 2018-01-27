function [finaldata,normdata]= Surf_Fit(Data1,Data2, data_MZ, data_RT, thr)

% load('feature_ranks.mat')


Data=[Data1, Data2];
[feature_ranks]=rank_with_ttest(Data1,Data2);
BiomarkersSel=feature_ranks(find(feature_ranks(:,2)<=27),1);
Housekeeping=setxor(1:size(Data,1),BiomarkersSel);
data_MZ_HK=data_MZ(Housekeeping,:);
data_RT_HK=data_RT(Housekeeping,:);

%% Scaling Data 

% % Median Scaling
% MedianData=median(Data);
% MedianData=repmat(MedianData,651,1);
% normdata=Data./MedianData;
% 
% 
% %MAD scaling
% MadData=mad(Data);
% MadData=repmat(MadData,651,1);
% normdata=Data./MadData;



%MinMax Scaling
% [normdata,minn,maxx]=minmax_Scaling(Data,0);
[normdata,standev,means]=Pareto_Scaling(Data,2);


%%  Dividing by Median
medians=[];
testmedian=[];
testmediann=[];
normdataa=normdata+max(median(normdata,2));
% normdataa=normdata;
medians=median(normdataa,2);
medians=abs(medians);
test=repmat(medians,1,size(Data,2));
testmediann=(normdataa)./(test);
testmedian=testmediann(Housekeeping,:);


%% Surface Fitting
parfor i =1:20
    [ff]=fittingbysample(testmedian(:,i),data_MZ_HK(:,i),data_RT_HK(:,i),thr);
    fff{i,1}=ff{1,1};

%      figure()
%     plot(ff, [data_MZ_HK(:,i),data_RT_HK(:,i)], testmedian(:,i))
end


[finaldataN]=estimateData_afterFitting(fff,data_MZ,data_RT);



%% get Final Data
finaldata=(normdataa(:,1:20))./abs(finaldataN);

