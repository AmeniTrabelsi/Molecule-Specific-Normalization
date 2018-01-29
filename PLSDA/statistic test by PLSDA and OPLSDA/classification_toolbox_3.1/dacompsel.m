function res = dacompsel(X,class,cv_type,cv_groups,class_prob,method,scal,max_comp)

% Selection of the optimal number of components for PCA coupled with DA
% by means of cross-validation
%
% res = dacompsel(x,y,cv_type,cv_groups,class_prob,method,scal,max_comp)
%
% input:
% X             dataset [samples x variables]
% class         class vector [samples x 1]
% cv_type       type of cross validation
%               'vene' for venetian blinds'
%               'cont' for contiguous blocks
% cv_groups     number of cv groups
%               if num_can == n (number of samples): leave-one-out
% class_prob    if prob = 1 equal probability, if prob = 2 proportional prob.
% method        'linear' or 'quadratic'
% scal          scaling method:
%               if scal = 'none' -> no scaling
%               if scal = 'cent' -> centering
%               if scal = 'auto' -> autoscaling
% max_comp      maximum number of components to be evaluated
%
% output:
% res is a structure array with fields:
% er            error rate in cross-validation for each component
% ner           not-error rate in cross-validation for each component
% not_ass:      ratio of not-assigned samples for each component
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
if r > max_comp
    r = max_comp;
end

hwait = waitbar(0,'cross validating models','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(hwait,'canceling',0)

for k = 1:r
    waitbar(k/r)
    out = dacv(X,class,cv_type,cv_groups,class_prob,method,k,scal);
    res.er(k) = out.class_param.er;
    res.ner(k) = out.class_param.ner;
    res.not_ass(k) = out.class_param.not_ass;
end

delete(hwait)