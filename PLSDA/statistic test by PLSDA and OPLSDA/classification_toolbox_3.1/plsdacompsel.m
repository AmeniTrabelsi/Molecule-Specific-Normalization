function res = plsdacompsel(X,class,pret_type,cv_type,cv_groups,assign_method)

% PLSDA selection of the optimal number of components 
% by means of cross-validation
%
% res = plsdacompsel(X,class,pret_type,cv_type,cv_groups,assign_method)
%
% input:
% X                 dataset [samples x variables]
% class             class vector [samples x 1]
% pret_type         scaling method:
%                   if scal = 'none' -> no scaling
%                   if scal = 'cent' -> centering
%                   if scal = 'auto' -> autoscaling
% cv_type           type of cross validation
%                   'vene' for venetian blinds'
%                   'cont' for contiguous blocks
% cv_groups         number of cv groups
%                   if num_can == n (number of samples): leave-one-out
% assign_method     assignation method
%                   'bayes' samples are assigned on thresholds based on Bayes Theorem
%                   'max' samples are assigned to the class with maximum yc
%
% output:
% res is a structure array with fields:
% er:           error rate in cross-validation for each component
% ner:          not-error rate in cross-validation for each component
% not_ass:      ratio of not-assigned samples for each component
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

[n,p] = size(X);
r = min(n,p);
if r > 20
    r = 20;
end
if r > 2
    r = r - 1;
end
hwait = waitbar(0,'cross validating models','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(hwait,'canceling',0)

for k = 1:r
    waitbar(k/r)
    out = plsdacv(X,class,k,pret_type,cv_type,cv_groups,assign_method);
    res.er(k) = out.class_param.er;
    res.ner(k) = out.class_param.ner;
    res.not_ass(k) = out.class_param.not_ass;
end

delete(hwait)
