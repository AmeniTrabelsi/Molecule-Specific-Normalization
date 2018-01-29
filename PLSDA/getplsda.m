function [plspara]= getplsda(data, gnd)


% addpath('./PLSDA/statistic test by PLSDA and OPLSDA/classification_toolbox_3.1');
% addpath('./PLSDA/statistic test by PLSDA and OPLSDA/libPLS_1.6');

%%% normalization
norm_data = [];
delete_p = [];
for i = 1:size(data,1)
    if all(data(i, :) == 0)
        delete_p = [delete_p; i];
    else
        a = mean(data(i,:));
        b = std(data(i,:));
        norm_data = [norm_data; (data(i,:)-a)/b];
    end
end
data(delete_p, :) = [];

A=20;K=5;method='center';
CV = plscv(data',gnd',A,K,method);  %%%gnd are the label information corresponding to the group information
ncomp = CV.optLV;
if ncomp <= 2
    ncomp = 3;
elseif ncomp > 20
    ncomp = 20;
end

plspara = plsdafit(norm_data', gnd' ,ncomp, 'scal', 'bayes', 1);
plspara.vip = getVIP(size(norm_data,1), ncomp, plspara.W, plspara.expvar(:,2)');

