function [model] = pca_project(Xnew,model);

% [model] = pca_project(X,num_comp,scal);
%
% pca_project is a module for PCA analysis.
% it projects new samples in the defined PCA score space
%
% ------------ INPUT ---------------------------------------------------
% Xnew       data matrix [n x p]  n objects, p variables
% model      defined PCA model the samples will be projected in
%
% ------------ OUTPUT --------------------------------------------------
% model, structure with fields:
% Tnew       score matrix [n x num_comp] of the projected samples
% Thotpred   T2 Hotelling [n x 1]
% Qrespred   Q residuals [n x 1]
%
% 
% Classification toolbox for MATLAB
% version 3.1 - October 2013
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% http://michem.disat.unimib.it/chm/

X_in = test_pretreatment(Xnew,model.set.param);
T = X_in*model.L;

% T2 hotelling
I = zeros(size(T,2),size(T,2));
for i=1:size(T,2)
    I(i,i) = model.E(i);
end
for i=1:size(T,1)
    Thot(i) = T(i,:)*inv(I)*T(i,:)';
end

% Q residuals
Xmod = T*model.L';
Err = X_in - Xmod;
for i=1:size(T,1)
    Qres(i) = Err(i,:)*Err(i,:)';
end

model.Tpred = T;
model.Thotpred = Thot;
model.Qrespred = Qres;
