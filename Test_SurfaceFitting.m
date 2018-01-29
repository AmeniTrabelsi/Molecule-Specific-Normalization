load('new_data_comp.mat')
load('data_MZ.mat')
load('data_RT.mat')


 Modified_data=new_data_comp;
 group0=[1:10];
 group5=[11:20];
 
% % If we add noise
% 
per=2;
noise1=(new_data_comp(:,3).*((per*rand(651,1))./10));
Modified_data(:,3)=new_data_comp(:,3)+ noise1;
noise2=(new_data_comp(:,6).*((per*rand(651,1))./10));
Modified_data(:,6)=new_data_comp(:,6)+ noise2;
% noise3=(new_data_comp(:,12).*((per*rand(651,1))./10));
% Modified_data(:,12)=new_data_comp(:,12)+ noise3;



%% Try Surface Fitting
%Boxplot threshold thr
thr=2;

[finaldata,~]=Surf_Fit(Modified_data(:,group0), Modified_data(:,group5), data_MZ, data_RT, thr);


%% Try Other Normalization Methods

[Contrast_g05] = Contrast(Modified_data)
[CyclicLoess_g05] = CyclicLoess(Modified_data)
[Quantile_g05] = Quantile(Modified_data)

%% Feature Selection Tests

% % Rank with t-test
% [ranks]= rank_with_ttest(Modified_data(:,group0), Modified_data(:,group5))
% [ranks2]= rank_with_ttest(new_data_comp(:,group0), new_data_comp(:,group5))
% [ranksf]= rank_with_ttest(finaldata(:,group0), finaldata(:,group5))
% [ranksContrast]= rank_with_ttest(Contrast_g05(:,group0), Contrast_g05(:,group5))
% [ranksCyclicLoess]= rank_with_ttest(CyclicLoess_g05(:,group0) ,CyclicLoess_g05(:,group5))
% [ranksQuantile]= rank_with_ttest(Quantile_g05(:,group0), Quantile_g05(:,group5))

%or with PLSDA

[ranks]= rank_with_plsda(Modified_data(:,group0), Modified_data(:,group5))
[ranks2]= rank_with_plsda(new_data_comp(:,group0), new_data_comp(:,group5))
[ranksf]= rank_with_plsda(finaldata(:,group0), finaldata(:,group5))
[ranksContrast]= rank_with_plsda(Contrast_g05(:,group0), Contrast_g05(:,group5))
[ranksCyclicLoess]= rank_with_plsda(CyclicLoess_g05(:,group0) ,CyclicLoess_g05(:,group5))
[ranksQuantile]= rank_with_plsda(Quantile_g05(:,group0), Quantile_g05(:,group5))


%% Plot ROC Curve of all the methods

figure();[treshP tp auc]=calc_roc_t(ranks,27,1,'b')
hold on;[treshP tp auc]=calc_roc_t(ranksf,27,1,'m')
hold on;[treshP tp auc]=calc_roc_t(ranksContrast,27,1,'c')
hold on;[treshP tp auc]=calc_roc_t(ranksCyclicLoess,27,1,'r')
hold on;[treshP tp auc]=calc_roc_t(ranksQuantile,27,1,'g')
title('ROC for different Normalization methods using PLSDA on Group 0 and 5')
legend('No Normalization','MSN','Contrast','Cyclic Loess','Quantile')
