function F=roccurve(ypred,yreal,flag)
%+++ yreal: with elements 1 or -1;
%+++ ypred: real values.
%+++ flag: 1: plot
%          0: no plot.
%+++ June 11,2008.

if nargin<3;flag=1;end;

yreal=sign(yreal);
thitamin=min(ypred);thitamax=max(ypred);
K=128;
if K>size(yreal,1);K=size(yreal,1);end

thita=linspace(thitamin-0.000001,thitamax+0.000001,K);
Result=zeros(K,2);
for i=1:K
  r=sesp(ypred-thita(i),yreal);
  Result(i,:)=[r.specificity r.sensitivity];    
end
auc=abs(trapz(Result(:,1),Result(:,2)));
if flag==1
  plot(1-Result(:,1),Result(:,2));
  xlabel('1-specificity');ylabel('sensitivity');
end
r=sesp(yreal,ypred);
%+++ OUTPUT
F.value=Result;
F.sensitivity=r.sensitivity;
F.specificity=r.specificity;
F.accuracy=r.accuracy;
F.AUC=auc; % area under curve.

