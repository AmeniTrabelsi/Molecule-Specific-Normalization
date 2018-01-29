function assigned_class = plsdafindclass(yc,class_thr);

% assign samples for PLSDA on the basis of thresholds and calculated responses
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

nobj = size(yc,1);
nclass = size(yc,2);
for i = 1:nobj
    pred = yc(i,:);
    chk_ass = zeros(1,nclass);
    for c = 1:nclass
        if pred(c) > class_thr(c); chk_ass(c) = 1; end;
    end
    if length(find(chk_ass)) == 1
        assigned_class(i) = find(chk_ass);
    else
        assigned_class(i) = 0;
    end
end