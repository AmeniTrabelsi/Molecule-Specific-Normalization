function pca_plot(model,comp,label_samples,label_var)

% pca_plot is a module for PCA analysis.
% it makes score and loading plots of the defined PCA model
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

L = model.L;
T = model.T;
Thot = model.Thot;
Qres = model.Qres;
if isfield(model,'Tnew') Tnew = model.Tnew; end

x = comp(1);
y = comp(2);
lab_x = (['PC ' num2str(x) ' - EV=' num2str(round(model.exp_var(x)*100)) '%']);
lab_y = (['PC ' num2str(y) ' - EV=' num2str(round(model.exp_var(y)*100)) '%']);

if nargin < 3
    label_samples = [];
    label_var = [];
end

figure
hold on
scatter(T(:,x),T(:,y),'.k')
if length(label_samples) == 0
    for j=1:size(T,1); text(T(j,x),T(j,y),num2str(j)); end
else
    for j=1:size(T,1); text(T(j,x),T(j,y),label_samples(j)); end
end
if isfield(model,'Tnew') scatter(Tnew(:,x),Tnew(:,y),'.r'); end
title('score plot')
xlabel(lab_x)
ylabel(lab_y)
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
if ylim(1) < 0 & ylim(2) > 0; line(xlim,[0 0],'Color','k','LineStyle',':'); end
if xlim(1) < 0 & xlim(2) > 0; line([0 0],ylim,'Color','k','LineStyle',':'); end
hold off

figure
hold on
scatter(L(:,x),L(:,y),'.r')
if length(label_var) == 0
    for j=1:size(L,1); text(L(j,x),L(j,y),num2str(j)); end
else
    for j=1:size(L,1); text(L(j,x),L(j,y),label_var(j)); end
end
title('loading plot')
xlabel(lab_x)
ylabel(lab_y)
xlabel(lab_x)
ylabel(lab_y)
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
if ylim(1) < 0 & ylim(2) > 0; line(xlim,[0 0],'Color','k','LineStyle',':'); end
if xlim(1) < 0 & xlim(2) > 0; line([0 0],ylim,'Color','k','LineStyle',':'); end
hold off

% Q and T statistic
figure
tlim = model.set.tlim;
qlim = model.set.qlim;
hold on
plot(Thot,Qres,'.k')
if length(label_samples) == 0
    for j=1:size(T,1); text(Thot(j),Qres(j),num2str(j)); end;
else
    for j=1:size(T,1); text(Thot(j),Qres(j),label_samples(j)); end;
end
ylim = get(gca, 'YLim');
line([tlim tlim],ylim,'Color','r','LineStyle',':')
xlim = get(gca, 'XLim');
line(xlim,[qlim qlim],'Color','r','LineStyle',':')
title(['T2 and Q statistic'])
xlabel('T2')
ylabel('Q')