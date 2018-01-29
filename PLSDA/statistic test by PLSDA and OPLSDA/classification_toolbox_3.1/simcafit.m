function model = simcafit(X,class,num_comp,pret_type,assign_method)

% fit SIMCA model
%
% model = simcafit(X,class,num_comp,pret_type,assign_method)
%
% input:
% X                 dataset [samples x variables]
% class             class vector [samples x 1]
% comp              number of components for each class model [1 x classes]
% pret_type         data pretreatment 
%                   'cent' cenering
%                   'scal' variance scaling
%                   'auto' for autoscaling (centering + variance scaling)
%                   'rang' range scaling (0-1)
%                   'none' no scaling
% assign_method     assignation method
%                   'prob', samples are assigned to the class with
%                   probability higher than 50%. Probabilities are
%                   calculated on the basis of Q resulduals and T2 Hotelling
%                   'dist', samples are always assigned to the closest
%                   class. Distances are calculated merging the reduced Q
%                   resulduals and T2 Hotelling
%
% output:
% model structure containing
% T                 Scores for each class model [samples x comp] 
% Thot              T2 Hotelling for each class model [n x 1]
% Thot_reduced      reduced T2 Hotelling for each class model [n x 1]
% tlim              T2 Hotelling confidence limit for each class model [classes x 1]
% Qres              Q residuals for each class model [n x 1]
% Qres_reduced      reduced Q residuals for each class model [n x 1]
% qlim              Q residuals confidence limit for each class model [classes x 1]
% modelpca          structure containing the class pca models {classes x 1}
% class_calc_dist   class assignations calculated on the closests class [n x 1]
% dist              distances based on Qres_reduced and Thot_reduced
% prob              class probabilities [n x classes]
% class_calc_prob   class assignations calculated on the class with prob higher than 50% [n x 1]
% class_calc        calculated class [samples x 1]
% class_param       structure with error rate, confusion matrix, specificity, sensitivity, precision
% set               structure with model settings
%
% The main routine is class_gui
%
% Reference for the calculation of SIMCA distances and probabilities:
% Reference: http://wiki.eigenvector.com/index.php?title=Sample_Classification_Predictions
%
% Note that a detailed HTML help is provided with the toolbox.
% See the HTML HELP files (help.htm) for futher details and examples
%
% Classification toolbox for MATLAB
% version 3.1 - October 2013
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% http://michem.disat.unimib.it/chm/

for g=1:max(class)
    train = X(find(class==g),:);
    test = X(find(class~=g),:);
    modelpca = pca_model(train,num_comp(g),pret_type,0);
    modelpca = pca_project(test,modelpca);
    % scores
    T = zeros(size(X,1),num_comp(g));
    T(find(class==g),:) = modelpca.T;
    T(find(class~=g),:) = modelpca.Tpred;
    model.T{g} = T;
    % t hotelling
    Thot = zeros(size(X,1),1);
    Thot(find(class==g),:) = modelpca.Thot';
    Thot(find(class~=g),:) = modelpca.Thotpred';
    model.Thot{g} = Thot;
    model.Thot_reduced{g} = Thot/modelpca.set.tlim;
    model.tlim(g) = modelpca.set.tlim;
    % q residuals
    Qres = zeros(size(X,1),1);
    Qres(find(class==g),:) = modelpca.Qres';
    Qres(find(class~=g),:) = modelpca.Qrespred';
    model.Qres{g} = Qres;
    model.Qres_reduced{g} = Qres./modelpca.set.qlim;
    model.qlim(g) = modelpca.set.qlim;
    % pca model
    model.modelpca{g} = modelpca;
end

% always assign samples
for i=1:size(X,1)
    for g=1:max(class)
        Q = model.Qres_reduced{g};
        T = model.Thot_reduced{g};
        q = Q(i);
        t = T(i);
        d(i,g) = (q^2 + t^2)^0.5;
    end
    [v,c] = min(d(i,:));
    model.class_calc_dist(i,1) = c;
end
model.dist = d;

% assign samples with probabilities
for g=1:max(class)
    % Thot prob
    Tclass = model.Thot{g};
    F = Tclass.*(size(X,1) - num_comp(g))/(num_comp(g)*(size(X,1)-1));
    prob_thot(:,g) = 1 - fcdf(F,num_comp(g),size(X,1) - num_comp(g));
    % Qres prob
    A = model.Qres{g};
    [T,Eg,L] = svd( model.modelpca{g}.Err,0);
    reseig = diag(Eg).^2/(length(find(class == g)) - 1);
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
    for g=1:max(class)
        p(i,g) = min([prob_thot(i,g) prob_qres(i,g)]);
        simcathreshold = 0.95;
        phere = 0.5*(p(i,g)/(1 - simcathreshold));
        if phere > 1; phere = 1; end
        if phere < 0; phere = 0; end
        model.prob(i,g) = phere;
    end
    [p,c] = sort(-model.prob(i,:));
    p = -p;
    if length(find(p > 0.5)) == 1
        model.class_calc_prob(i,1) = c(1);
    else
        model.class_calc_prob(i,1) = 0;
    end
end
if strcmp(assign_method,'dist')
    model.class_calc = model.class_calc_dist;
else
    model.class_calc = model.class_calc_prob;
end
model.class_param = calc_class_param(model.class_calc,class);
model.set.class_true = class;
model.set.scal = pret_type;
model.set.assign_method = assign_method;