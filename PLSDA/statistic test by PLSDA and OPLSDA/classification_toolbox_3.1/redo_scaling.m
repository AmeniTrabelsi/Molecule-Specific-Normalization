function X = redo_scaling(X_scal,param)

% redo scaling (from scaled data to original data)
%
% input:
% X_scal:   pretreated data matrix [samples x variables]
% param:    output data structure from data_pretreatment routine
%
% output:
% X:        data matrix [samples x variables]
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

a = param.a;
s = param.s;
m = param.m;
M = param.M;
pret_type = param.pret_type;

if strcmp(pret_type,'cent')
    for i=1:size(X_scal,1)
        for j=1:size(X_scal,2)
            X(i,j) = X_scal(i,j) + a(j);
        end
    end
elseif strcmp(pret_type,'scal')
    for i=1:size(X_scal,1)
        for j=1:size(X_scal,2)
            X(i,j) = X_scal(i,j)/s(j);
        end
    end
elseif strcmp(pret_type,'auto')
    for i=1:size(X_scal,1)
        for j=1:size(X_scal,2)
            X(i,j) = X_scal(i,j)*s(j) + a(j);
        end
    end
elseif strcmp(pret_type,'rang')
    for i=1:size(X_scal,1)
        for j=1:size(X_scal,2)
            X(i,j) = X_scal(i,j)*(M(j) - m(j)) + m(j);
        end
    end
else
    X = X_scal;
end
