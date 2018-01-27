function [treshP tp auc]=calc_roc_t(rank,numBio,tresh,s)
%rank contains two collumns, 1-feature numbers and 2-rank of features
numFea=size(rank,1);
%figure
hold on
xlabel('FPR')
ylabel('TPR')
    FR=sortrows(rank,2);
    x(1,1:numFea)=1:numFea;
    for i=1:numFea
        x(2,i)=sum(FR(1:i,1)<=numBio)/numBio;   
        x(3,i)=sum(FR(1:i,1)>numBio)/(numFea-numBio);
    end
    axis([0 0.1 0 1]);
    plot(x(3,:),x(2,:),s,'LineWidth',2)
    %title(['AUC',num2str( trapz(x(3,:),x(2,:)))]);
    %trapz(x(3,:),x(2,:))
    %hold off
    treshP=min(find(x(3,:)>=tresh));
    auc=trapz(x(3,1:treshP),x(2,1:treshP));
    tp=x(2,min(find(x(3,:)>tresh)));
end
