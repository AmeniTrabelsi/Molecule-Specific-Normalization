function pred = plsdapred(X,model)

% apply PLSDA model on unknow samples
%
% pred = plsdapred(X,model)
%
% input:
% X                 dataset [samples x variables]
% model             plsda model (calculated with plsdafit routine)
%
% output:
% pred is a structure containing:
% yc                predicted response [samples x classes]
% class_pred        predicted class vector [samples x 1]
% H                 leverages [samples x 1]
% T                 predicted scores [samples x comp]
% Thot              T2 Hotelling [samples x 1]
% Qres              Q residuals [samples x 1]
%
% based on Frans Van Den Berg mypls routine
% http://www.models.kvl.dk/
%
% The main routine is class_gui
%
% Note that a detailed HTML help is provided with the toolbox.
% See the HTML HELP files (help.htm) for futher details and examples
%
% Classification toolbox for MATLAB
% version 3.1 - October 2013
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% http://michem.disat.unimib.it/chm/

x = X;
W = model.W;
Q = model.Q;
P = model.P;
nF = size(model.T,2);
T = model.T;
xscal = test_pretreatment(x,model.set.px);

yscal_c = 0;
for k = 1:nF
    Ttest(:,k) = xscal*W(:,k)/(W(:,k)'*W(:,k));
    yscal_c = yscal_c + Ttest(:,k)*Q(:,k)';
    xscal = xscal - Ttest(:,k)*P(:,k)';
end

yc = redo_scaling(yscal_c,model.set.py);

if strcmp(model.set.assign_method,'max')
    [non,assigned_class] = max(yc');
else
    assigned_class = plsdafindclass(yc,model.set.thr);
end
pred.class_pred = assigned_class';
pred.yc = yc;

% leverages
xscal = test_pretreatment(x,model.set.px);
pred.H = diag(Ttest*pinv(T'*T)*Ttest');
pred.T = Ttest;

% T hot
fvar = sqrt(1./(diag(T'*T)/(size(T,1) - 1)));
pred.Thot = sum((Ttest*diag(fvar)).^2,2);

% Qres
Xmod = Ttest*P';
E = xscal - Xmod;
for i=1:size(X,1)
    pred.Qres(i) = E(i,:)*E(i,:)';
end