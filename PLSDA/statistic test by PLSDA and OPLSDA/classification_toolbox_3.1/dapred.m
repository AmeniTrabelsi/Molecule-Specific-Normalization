function pred = dapred(X,model)

% prediction of new samples with Discriminant Analysis
%
% pred = dapred(X,model)
%
% input:
% X             dataset [samples x variables]
% model         Discriminant Analysis model
%
% output:
% pred is a structure conyaining
% class_pred    predicted class vector [samples x 1]
% S             scores on canonical variables [samples x G-1], only for lda
% T             scores on PCA if da is calculate on principal components
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

class_prob = model.settings.class_prob;
method = model.settings.method;
Xtrain = model.settings.Xtrain;
num_comp = model.settings.num_comp;
class_train = model.settings.class_true;
nclass = max(class_train);

if class_prob == 2
    for g = 1:nclass 
        obj_cla(g)  = length(find(class_train == g));
    end
    prob = obj_cla/size(Xtrain,1);
end

if num_comp > 0
    modelpca = pca_project(X,model.settings.modelpca);
    Xin = modelpca.Tpred;
else
    Xin = X;
end

% if linear and not with PCs check for pooled estimate of covariance
doit = 1;
if strcmp('linear',method) & num_comp == 0
    doit = pec(Xtrain,class_train);
end

if doit
    % prediction
    if class_prob == 1
        [class_pred] = classify(Xin,Xtrain,class_train,method);
    else
        [class_pred] = classify(Xin,Xtrain,class_train,method,prob);
    end
    % prediction of scores on canonical variables only for lda
    if strcmp('linear',method)
        [a,param] = data_pretreatment(Xtrain,'cent');
        Xin_cent = test_pretreatment(Xin,param);
        pred.S = Xin_cent*model.L;
    end
else
    class_pred = ones(size(Xin,1),1);
end

pred.class_pred = class_pred;
if num_comp > 0
    pred.T = modelpca.Tpred;
end

% -------------------------------------------------------------------------
function doit = pec(X,class)

for g = 1:max(class)
    gmeans(g,:) = mean(X(find(class == g),:));
end
% Pooled estimate of covariance
[Q,R] = qr(X - gmeans(class,:), 0);
R = R / sqrt(size(X,1) - max(class)); % SigmaHat = R'*R
s = svd(R);
if any(s <= eps^(3/4)*max(s))
    doit = 0;
    disp('The pooled covariance matrix of training samples must be positive definite. model not calculated');
else
    doit = 1;
end