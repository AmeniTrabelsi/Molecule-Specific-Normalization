function [model] = pca_model(X,num_comp,scal,do_plot);

% [model] = pca_model(X,num_comp,scal);
%
% pca_model is a module for PCA analysis.
% it calculates the PCA model on a data matrix
%
% ------------ INPUT ---------------------------------------------------
% X          data matrix [n x p]  n objects, p variables
% num_comp   number of principal components you want to retain
% scal       scaling method:
%            if scal = 'none' -> no scaling
%            if scal = 'cent' -> centering
%            if scal = 'auto' -> autoscaling
% do_plot    optional, if plot = 1 does plots
%
% ------------ OUTPUT --------------------------------------------------
% model, structure with fields:
% exp_var    explained variance [num_comp x 1]
% cum_var    cumulative explained variance [num_comp x 1]
% T          score matrix [n x num_comp]
% L          loading matrix [p x num_comp]
% E          eigenvalues [num_comp x 1]
% Thot       T2 Hotelling [n x 1]
% Qres       Q residuals [n x 1]
% set        structure array with num_comp, scal, mean and st. dev, raw data
%            limits of q and t statistics
%
% Classification toolbox for MATLAB
% version 3.1 - October 2013
% Davide Ballabio
% Milano Chemometrics and QSAR Research Group
% http://michem.disat.unimib.it/chm/

if nargin < 4
    do_plot = 0;
end

[n,p] = size(X);
ran = min(size(X,1),size(X,2));
if num_comp > ran
    num_comp = ran;
end

[X_in,param] = data_pretreatment(X,scal);

[T,E,L] = svd(X_in,0);     % diagonalisation
eigmat = E;
Tmat = T;
E = diag(E).^2/(n-1);      % eigenvalues
exp_var = E/sum(E);
E = E(1:num_comp);
exp_var = exp_var(1:num_comp);
for k=1:num_comp; cum_var(k) = sum(exp_var(1:k)); end;

L = L(:,1:num_comp);       % loadings and scores
T = X_in*L;
T = T(:,1:num_comp);

% T2 hotelling
I = zeros(size(T,2),size(T,2));
for i=1:size(T,2)
    I(i,i) = E(i);
end
for i=1:size(T,1)
    Thot(i) = T(i,:)*inv(I)*T(i,:)';
end

% Q residuals
Xmod = T*L';
Err = X_in - Xmod;
for i=1:size(T,1)
    Qres(i) = Err(i,:)*Err(i,:)';
end

% T2 and Q limits
[tlim,qlim] = calc_qt_pca(num_comp,size(X,1),Err);

% save results
model.exp_var = exp_var;
model.cum_var = cum_var';
model.E = E;
model.L = L;
model.T = T;
model.Err = Err;
model.eigmat = eigmat;
model.Tmat = Tmat;
model.Qres = Qres;
model.Thot = Thot;
model.set.tlim = tlim;
model.set.qlim = qlim;
model.set.num_comp = num_comp;
model.set.param = param;

% plots
if(do_plot)
    figure
    marker_size = 15;
    subplot(3,1,1)
    hold on
    plot(E,'k')
    plot(E,'.k','MarkerSize',marker_size)
    ylim = get(gca, 'YLim');
    axis([0.6 (num_comp + 0.4) ylim(1) ylim(2)])
    set(gca,'xtick',[1:num_comp]);
    hold off
    ylabel('eigenvalues')
    
    subplot(3,1,2)
    hold on
    plot(exp_var*100,'k')
    plot(exp_var*100,'.k','MarkerSize',marker_size)
    ylim = get(gca, 'YLim');
    axis([0.6 (num_comp + 0.4) ylim(1) ylim(2)])
    set(gca,'xtick',[1:num_comp]);
    hold off
    ylabel('explained variance (%)')
    
    subplot(3,1,3)
    hold on
    plot(cum_var*100,'k')
    plot(cum_var*100,'.k','MarkerSize',marker_size)
    axis([0.6 (num_comp + 0.4) 0 100])
    set(gca,'xtick',[1:num_comp]);
    hold off
    ylabel('cumulative variance (%)')
    xlabel('principal components')
end