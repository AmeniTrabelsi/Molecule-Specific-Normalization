function cv = plsdacv(X,class,comp,pret_type,cv_type,cv_groups,assign_method)

% cross-validation for PLSDA
%
% cv = plsdacv(X,class,comp,pret_type,cv_type,cv_groups,assign_method)
%
% X                 dataset [samples x variables]
% class             class vector [samples x 1]
% comp              number of components
% pret_type         scaling method:
%                   if scal = 'none' -> no scaling
%                   if scal = 'cent' -> centering
%                   if scal = 'auto' -> autoscaling
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
%                   'bayes' samples are assigned on thresholds based on Bayes Theorem
%                   'max' samples are assigned to the class with maximum yc
%
% output:
% cv structure containing
% class_pred    predicted class vector [samples x 1] in cross-validation
% rmsec         root mean squared error inc ross-validation (1 x g) 
% class_param   contains error rate, confusion matrix, specificity, sensitivity, precision
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

y = class;
x = X;
if size(y,2) == 1
    tmp_unfold = zeros(size(x,1),max(y));
    for g=1:max(y)
        tmp_unfold(find(y==g),g) = 1;
    end
    y = tmp_unfold;
end
nobj=size(x,1);
if strcmp(cv_type,'boot')
    out_bootstrap = zeros(nobj,1);
    assigned_class = [];
    class_true = [];
    for i=1:cv_groups
        out = ones(nobj,1);
        whos_in = [];
        for k=1:nobj
            r = ceil(rand*nobj);
            whos_in(k) = r;
        end
        out(whos_in) = 0;
        % counters for left out samples
        boot_how_many_out(i)=length(find(out == 1));
        out_bootstrap(find(out == 1)) = out_bootstrap(find(out == 1)) + 1;
        
        x_out = x(find(out == 1),:);
        y_out = y(find(out == 1));
        x_in  = x(whos_in,:);
        y_in  = y(whos_in,:);
        
        model = plsdafit(x_in,y_in,comp,pret_type,assign_method,0);
        pred = plsdapred(x_out,model);
        assigned_class = [assigned_class; pred.class_pred];
        class_true = [class_true; class(find(out == 1))];
        for g=1:size(y,2)
            rmsec(g) = NaN;
            quantitative_class(g) = NaN;
        end
    end
    class = class_true;
    assigned_class = assigned_class';
elseif strcmp(cv_type,'rand')
    assigned_class = [];
    out_rand = zeros(nobj,1);
    perc_in = 0.8;
    take_in = round(nobj*perc_in);
    class_true = [];
    for i=1:cv_groups
        out = ones(nobj,1);
        whos_in = randperm(nobj);
        whos_in = whos_in(1:take_in);
        out(whos_in) = 0;
        % counters for left out samples
        out_rand(find(out == 1)) = out_rand(find(out == 1)) + 1;
        
        x_out = x(find(out == 1),:);
        y_out = y(find(out == 1));
        x_in  = x(whos_in,:);
        y_in  = y(whos_in,:);
        
        model = plsdafit(x_in,y_in,comp,pret_type,assign_method);
        pred = plsdapred(x_out,model);
        assigned_class = [assigned_class; pred.class_pred];
        class_true = [class_true; class(find(out == 1))];
        for g=1:size(y,2)
            rmsec(g) = NaN;
            quantitative_class(g) = NaN;
        end
    end
    class = class_true;
    assigned_class = assigned_class';
else
    quantitative_class = zeros(nobj,size(y,2));
    class_pred = zeros(nobj,1);
    obj_in_block = fix(nobj/cv_groups);
    left_over = mod(nobj,cv_groups);
    st = 1;
    en = obj_in_block;
    for i = 1:cv_groups
        in = ones(size(x,1),1);
        if strcmp(cv_type,'vene') % venetian blinds
            out = [i:cv_groups:nobj];
        else % contiguous blocks
            if left_over == 0
                out = [st:en];
                st =  st + obj_in_block;  en = en + obj_in_block;
            else
                if i < cv_groups - left_over
                    out = [st:en];
                    st =  st + obj_in_block;  en = en + obj_in_block;
                elseif i < cv_groups
                    out = [st:en + 1];
                    st =  st + obj_in_block + 1;  en = en + obj_in_block + 1;
                else
                    out = [st:nobj];
                end
            end
        end
        in(out) = 0;
        x_in = x(find(in),:);
        y_in = y(find(in),:);
        x_out = x(find(in == 0),:);
        model = plsdafit(x_in,y_in,comp,pret_type,assign_method,0);
        out = plsdapred(x_out,model);
        assigned_class(find(in == 0)) = out.class_pred;
        quantitative_class(find(in == 0),:) = out.yc;
    end
    for g=1:size(y,2)
        C = calc_reg_param(y(:,g),quantitative_class(:,g));
        rmsec(g) = C.RMSEC;
    end
end

class_param = calc_class_param(assigned_class',class);

cv.class_pred = assigned_class';
cv.class_param = class_param;
cv.yc = quantitative_class;
cv.rmsec = rmsec;
cv.settings.cv_groups = cv_groups;
cv.settings.cv_type = cv_type;
cv.settings.num_comp = comp;
cv.settings.scal = pret_type;