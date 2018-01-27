load('new_data_comp2.mat')

% If we add noise

 Modified_data=new_data_comp;
noise1=(new_data_comp(:,3).*((2*rand(651,1))./10));
 %noise1=((1*max(new_data_comp(:,3))*rand(1,1))./10);
Modified_data(:,3)=new_data_comp(:,3)+ noise1;
noise2=(new_data_comp(:,6).*((2*rand(651,1))./10));
Modified_data(:,6)=new_data_comp(:,6)+ noise2;

Modified_data=new_data_comp;
sel=randperm(651,floor(0.8*651));
Modified_data(sel,3)=1.3*new_data_comp(sel,3);
sel2=randperm(651,floor(0.8*651));
Modified_data(sel2,5)=1.3*new_data_comp(sel2,5);
sel3=randperm(651,floor(0.2*651));
Modified_data(sel3,7)=1.3*new_data_comp(sel3,7);


Modified_data=new_data_comp;
noise1=(new_data_comp(:,3).*((4*normrnd(651,1))./10));
 %noise1=((1*max(new_data_comp(:,3))*rand(1,1))./10);
Modified_data(:,3)=new_data_comp(:,3)+ noise1;
%% Try Surface Fitting

load('data_MZ.mat')
load('data_RT.mat')

[finaldata,normdata]=Surf_Fit(Modified_data(:,1:10), Modified_data(:,11:20), data_MZ, data_RT, 2);

[Yfcl Zfcl]=smoothh(Modified_data);

[Yfcl Zfcl]=smoothh(finaldata);
finaldataa(:,1:20)=[Yfcl Zfcl];

%% Try Other Normalization Methods

[Contrast_grp0_20perN] = Contrast(Modified_data)
[CyclicLoess_grp0_20perN] = CyclicLoess(Modified_data)
[Quantile_grp0_20perN] = Quantile(Modified_data)

%% Fusion Tests

[feature_ranks2 results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[finaldata(:,1:10) finaldata(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)
[feature_rankss results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[normdata(:,1:10) normdata(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)
[feature_ranksss results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[finaldataa(:,1:10) finaldataa(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)

[feature_ranksTruth results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[new_data_comp(:,1:10) new_data_comp(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)
[feature_ranksTruth results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[Modified_data(:,1:10) Modified_data(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)

[feature_ranksContrast results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[Contrast_grp0_20perN(:,1:10) Contrast_grp0_20perN(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)
[feature_ranksCyclicLoess results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[CyclicLoess_grp0_20perN(:,1:10) CyclicLoess_grp0_20perN(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)
[feature_ranksQuantile results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[Quantile_grp0_20perN(:,1:10) Quantile_grp0_20perN(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)
[feature_ranksTrimedCMean results truely_selected_biomarkers nts]=Ensemble_Feature_Selection([[TrimedCMean_grp0_10perN(:,1:10) TrimedCMean_grp0_10perN(:,11:20)]; zeros(1,10) 5*ones(1,10)]',27,11,0.01,10,0,5)


[ranks]= rank_with_ttest(Modified_data(:,1:10), Modified_data(:,11:20))
[ranks2]= rank_with_ttest(new_data_comp(:,1:10), new_data_comp(:,11:20))
[ranksf]= rank_with_ttest(finaldata(:,1:10), finaldata(:,11:20))
[ranksContrast]= rank_with_ttest(Contrast_grp0_20perN(:,1:10), Contrast_grp0_20perN(:,11:20))
[ranksCyclicLoess]= rank_with_ttest(CyclicLoess_grp0_20perN(:,1:10) ,CyclicLoess_grp0_20perN(:,11:20))
[ranksQuantile]= rank_with_ttest(Quantile_grp0_20perN(:,1:10), Quantile_grp0_20perN(:,11:20))

[ranks]= rank_with_plsda(Modified_data(:,1:10), Modified_data(:,11:20))
[ranks2]= rank_with_plsda(new_data_comp(:,1:10), new_data_comp(:,11:20))
[ranksf]= rank_with_plsda(finaldata(:,1:10), finaldata(:,11:20))
[ranksContrast]= rank_with_plsda(Contrast_grp0_20perN(:,1:10), Contrast_grp0_20perN(:,11:20))
[ranksCyclicLoess]= rank_with_plsda(CyclicLoess_grp0_20perN(:,1:10) ,CyclicLoess_grp0_20perN(:,11:20))
[ranksQuantile]= rank_with_plsda(Quantile_grp0_20perN(:,1:10), Quantile_grp0_20perN(:,11:20))

figure();[treshP tp auc]=calc_roc_t(ranks2,27,1,'b')
hold on;[treshP tp auc]=calc_roc_t(ranks,27,1,'k')
hold on;[treshP tp auc]=calc_roc_t(ranksf,27,1,'m')
hold on;[treshP tp auc]=calc_roc_t(ranksContrast,27,1,'c')
hold on;[treshP tp auc]=calc_roc_t(ranksCyclicLoess,27,1,'r')
hold on;[treshP tp auc]=calc_roc_t(ranksQuantile,27,1,'g')
%% Plot ROCs

figure();[treshP tp auc]=calc_roc_t(feature_ranks2,27,1,'b')
hold on;[treshP tp auc]=calc_roc_t(feature_rankss,27,1,'m')
hold on;[treshP tp auc]=calc_roc_t(feature_ranksss,27,1,'g')
hold on;[treshP tp auc]=calc_roc_t(feature_ranksTruth,27,1,'k')
hold on;[treshP tp auc]=calc_roc_t(feature_ranksContrast,27,1,'c')
hold on;[treshP tp auc]=calc_roc_t(feature_ranksCyclicLoess,27,1,'r')
hold on;[treshP tp auc]=calc_roc_t(feature_ranksQuantile,27,1,'g')

figure();
[tpr,fpr]=roc([ones(27,1); zeros(651-27,1)]',651-feature_ranks2(:,2)')
plot(fpr,tpr,'b')
hold on;
[tpr2,fpr2]=roc([ones(27,1); zeros(651-27,1)]',651-feature_ranksss(:,2)')
plot(fpr2,tpr2,'r')


%% Plot Biomarkers
figure()
for i=1 :25
        subplot(5,5,i)       % add first plot in 2 x 1 grid
        plot(1:10,finaldata(BiomarkersSel(i,1),1:10)-max(median(normdata,2)),'b',11:20,finaldata(BiomarkersSel(i,1),11:20)-max(median(normdata,2)),'b');
    hold on; plot(1:10,normdata(BiomarkersSel(i,1),1:10),'r',11:20,normdata(BiomarkersSel(i,1),11:20),'r');
%     hold on; plot(1:10,normdatatruth(BiomarkersSel(i,1),1:10),'g',11:20,normdatatruth(BiomarkersSel(i,1),11:20),'g');
%     hold on; plot(1:10,Quantile_grp0_20perN(BiomarkersSel(i,1),1:10),'g',11:20,Quantile_grp0_20perN(BiomarkersSel(i,1),11:20),'g');
    hold off;
        title(['Biomarker ' num2str(BiomarkersSel(i,1))])
end


%% Plot Normdata in 2D (MZ,data)
[~,order]=sort(data_MZ(:,1));
figure();plot(data_MZ(order,1),finaldataN(order,1));
hold on;plot(data_MZ(order,2),finaldataN(order,2));
hold on;plot(data_MZ(order,4),finaldataN(order,4));
% hold on;plot(data_MZ(order,6),finaldataN(order,6),'go');
hold on;plot(data_MZ(order,7),finaldataN(order,7));
hold on;plot(data_MZ(order,8),finaldataN(order,8));
hold on ;plot(data_MZ(order,3),finaldataN(order,3),'r*');
hold on ;plot(data_MZ(order,5),finaldataN(order,5));

[~,orderr]=sort(data_RT(:,1));
figure();plot(data_RT(orderr,1),finaldataN(orderr,1));
hold on;plot(data_RT(orderr,2),finaldataN(orderr,2));
hold on;plot(data_RT(orderr,4),finaldataN(orderr,4));
% hold on;plot(data_RT(orderr,6),finaldataN(orderr,6),'go');
hold on;plot(data_RT(orderr,7),finaldataN(orderr,7));
hold on;plot(data_RT(orderr,8),finaldataN(orderr,8));
hold on ;plot(data_RT(orderr,3),finaldataN(orderr,3),'r*');
hold on ;plot(data_RT(orderr,5),finaldataN(orderr,5));