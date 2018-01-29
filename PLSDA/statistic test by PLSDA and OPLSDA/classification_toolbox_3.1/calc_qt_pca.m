function [tlim,qlim] = calc_qt(comp,nobj,E);

% calc Q and T2 limits of confidence
% based on the presence of statistics_toolbox
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

% T2 limit
lev_conf = 0.95;
if license('test','statistics_toolbox')
    F = finv(lev_conf,comp,nobj-comp);
    tlim = comp*(nobj - 1)/(nobj - comp)*F;
else
    tlim = NaN;
end

% Q limit
[T,Eg,L] = svd(E,0);     % diagonalisation
eigmat = Eg;

reseig = diag(eigmat).^2/(nobj - 1);
m = length(reseig);
cl = 2*lev_conf*100-100;
t1 = sum(reseig(1:m,1));
t2 = sum(reseig(1:m,1).^2);
t3 = sum(reseig(1:m,1).^3);
if t1==0 
    qlim = 0;
else
    ca = sqrt(2)*erfinv(cl/100);
    h0 = 1-2*t1*t3/(3*(t2.^2)); 
    if h0<0.001; h0 = 0.001; end
    h1 = ca*sqrt(2*t2*h0.^2)/t1; h2 = t2*h0*(h0-1)/(t1.^2);
    qlim = t1*(1+h1+h2).^(1/h0);
end