function pred = simcapred(X,model)

% prediction with SIMCA model
%
% pred = simcapred(X,model)
%
% input:
% X                 dataset [samples x variables]
% model             structure containing the simca model calculated with
%                   simcafit routine
%
% output:
% model structure containing
% T                 Scores for each class model [samples x comp] 
% Thot              T2 Hotelling for each class model [n x 1]
% Thot_reduced      reduced T2 Hotelling for each class model [n x 1]
% Qres              Q residuals for each class model [n x 1]
% Qres_reduced      reduced Q residuals for each class model [n x 1]
% class_pred_dist   class assignations calculated on the closests class [n x 1]
% prob              class probabilities [n x classes]
% class_pred_prob   class assignations calculated on the class with prob higher than 50% [n x 1]
% class_pred        calculated class [samples x 1]
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

for g=1:max(model.set.class_true)
    modelpca = model.modelpca{g};
    modelpca = pca_project(X,modelpca);
    % scores
    pred.T{g} = modelpca.Tpred;
    % t hotelling
    pred.Thot{g} = modelpca.Thotpred';
    pred.Thot_reduced{g} = pred.Thot{g}./modelpca.set.tlim;
    % q residuals
    pred.Qres{g} = modelpca.Qrespred';
    pred.Qres_reduced{g} = pred.Qres{g}./modelpca.set.qlim;
end

% always assign samples
for i=1:size(X,1)
    for g=1:max(model.set.class_true)
        Q = pred.Qres_reduced{g};
        T = pred.Thot_reduced{g};
        q = Q(i);
        t = T(i);
        d(i,g) = (q^2 + t^2)^0.5;
    end
    [v,c] = min(d(i,:));
    pred.class_pred_dist(i,1) = c;
end

% assign samples with probabilities
for g=1:max(model.set.class_true)
    % Thot prob
    Tclass = pred.Thot{g};
    ntrain = length(model.class_calc);
    F = Tclass.*(ntrain - size(model.T{g},2))/(size(model.T{g},2)*(ntrain - 1));
    prob_thot(:,g) = 1 - fcdf(F,size(model.T{g},2),ntrain - size(model.T{g},2));
    % Qres prob
    A = pred.Qres{g};
    [T,Eg,L] = svd(model.modelpca{g}.Err,0);
    reseig = diag(Eg).^2/(length(find(model.set.class_true == g)) - 1);
    msize = length(reseig);
    t1 = sum(reseig(1:msize,1));
    t2 = sum(reseig(1:msize,1).^2);
    t3 = sum(reseig(1:msize,1).^3);
    h0 = 1-2*t1*t3/(3*(t2.^2));
    if h0<0.001; h0 = 0.001; end
    h1 = (A/t1).^h0;
    h2 = t2*h0*(h0-1)/t1^2;
    d = t1*(h1-1-h2)/(sqrt(2*t2)*h0);
    prob_qres(:,g) = 1 - 0.5*(1+erf(d/sqrt(2)));
end
for i=1:size(X,1)
    for g=1:max(model.set.class_true)
        p(i,g) = min([prob_thot(i,g) prob_qres(i,g)]);
        simcathreshold = 0.95;
        phere = 0.5*(p(i,g)/(1 - simcathreshold));
        if phere > 1; phere = 1; end
        if phere < 0; phere = 0; end
        pred.prob(i,g) = phere;
    end
    [p,c] = sort(-pred.prob(i,:));
    p = -p;
    if length(find(p > 0.5)) == 1
        pred.class_pred_prob(i,1) = c(1);
    else
        pred.class_pred_prob(i,1) = 0;
    end
end

if strcmp(model.set.assign_method,'dist')
    pred.class_pred = pred.class_pred_dist;
else
    pred.class_pred = pred.class_pred_prob;
end