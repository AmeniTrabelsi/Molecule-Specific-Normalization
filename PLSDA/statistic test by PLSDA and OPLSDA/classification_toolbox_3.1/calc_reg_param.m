function [res,RSS] = calc_reg_param(y_real,y_calc)

% calculation of regression parameters
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

n = length(y_real);

% TSS, Total Sum of Squares
TSS = sum((y_real - mean(y_real)).^2);

% RSS, Residual Sum of Squares
RSS = sum((y_real - y_calc).^2);

% R2, percentage of explained variance
R2 = 1 - RSS/TSS;

% R, coefficient of multiple correlation
R = R2^0.5;

% RMSEC, Root Mean Squared Error of Calibration
% also called SDEC, Standard Deviation Error in Calibration
RMSEC = (RSS/n)^0.5;

res.R2 = R2;
res.R = R;
res.RMSEC = RMSEC;