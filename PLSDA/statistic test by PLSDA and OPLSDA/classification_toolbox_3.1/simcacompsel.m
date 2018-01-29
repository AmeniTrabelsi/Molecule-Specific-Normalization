function res = simcacompsel(X,class,pret_type,cv_type,cv_groups,assign_method)

% SIMCA selection of the optimal number of components 
% by means of cross-validation
%
% res = simcacompsel(X,class,pret_type,cv_type,cv_groups,assign_method)
%
% input:
% X                 dataset [samples x variables]
% class             class vector [samples x 1]
% pret_type         data pretreatment 
%                   'cent' cenering
%                   'scal' variance scaling
%                   'auto' for autoscaling (centering + variance scaling)
%                   'rang' range scaling (0-1)
%                   'none' no scaling
% cv_type           type of cross validation
%                   'vene' for venetian blinds'
%                   'cont' for contiguous blocks
%                   'boot' for bootstrap with resampling
%                   'rand' for random sampling of 20% of samples
% cv_groups         number of cv groups
%                   if num_can == n (number of samples): leave-one-out
%                   if 'boot' or 'rand' are selected as cv_type, cv_groups 
%                   sets the number of iterations
% assign_method     assignation method
%                   'prob', samples are assigned to the class with
%                   probability higher than 50%. Probabilities are
%                   calculated on the basis of Q resulduals and T2 Hotelling
%                   'dist', samples are always assigned to the closest
%                   class. Distances are calculated merging the reduced Q
%                   resulduals and T2 Hotelling
%
% output:
% res is a structure array with fields:
% er:           error rate in cross-validation [classes x components]
% ner:          not-error rate in cross-validation [classes x components]
% not_ass:      ratio of not-assigned samples [classes x components]
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

[n,p] = size(X);
for g=1:max(class); c(g) = length(find(class==g));end
n = min(c);
r = min(n,p);
if r > 20
    r = 20;
end
if r > 2
    r = r - 1;
end
for g =1:max(class)
    hwait = waitbar(0,['cross validating models for class ' num2str(g)],'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(hwait,'canceling',0)
    for k = 1:r
        waitbar(k/r)
        num_comp = ones(1,max(class));
        num_comp(g) = k;
        out = simcacv(X,class,num_comp,pret_type,cv_type,cv_groups,assign_method);
        res.er(g,k) = out.class_param.er;
        res.ner(g,k) = out.class_param.ner;
        res.not_ass(g,k) = out.class_param.not_ass;
    end
    delete(hwait)
end
