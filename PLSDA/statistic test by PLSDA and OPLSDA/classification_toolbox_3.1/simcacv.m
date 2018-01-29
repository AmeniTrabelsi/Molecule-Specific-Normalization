function cv = simcacv(X,class,comp,pret_type,cv_type,cv_groups,assign_method)

% cross-validation for SIMCA
%
% cv = simcacv(X,class,comp,pret_type,cv_type,cv_groups,assign_method)
%
% X                 dataset [samples x variables]
% class             class vector [samples x 1]
% comp              number of components for each class model [1 x classes]
% pret_type         data pretreatment 
%                   'cent' cenering
%                   'scal' variance scaling
%                   'auto' for autoscaling (centering + variance scaling)
%                   'rang' range scaling (0-1)
%                   'none' no scaling
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
%                   'prob', samples are assigned to the class with
%                   probability higher than 50%. Probabilities are
%                   calculated on the basis of Q resulduals and T2 Hotelling
%                   'dist', samples are always assigned to the closest
%                   class. Distances are calculated merging the reduced Q
%                   resulduals and T2 Hotelling
%
% output:
% cv structure containing
% class_pred    predicted class vector [samples x 1] in cross-validation
% class_param   contains error rate, confusion matrix, specificity, sensitivity, precision
%
% The main routine is class_gui
%
% Reference for the calculation of SIMCA distances and probabilities:
% Reference: http://wiki.eigenvector.com/index.php?title=Sample_Classification_Predictions
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
        
        model = simcafit(x_in,y_in,comp,pret_type,assign_method);
        pred = simcapred(x_out,model);
        assigned_class = [assigned_class; pred.class_pred];
        class_true = [class_true; class(find(out == 1))];
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
        
        model = simcafit(x_in,y_in,comp,pret_type,assign_method);
        pred = simcapred(x_out,model);
        assigned_class = [assigned_class; pred.class_pred];
        class_true = [class_true; class(find(out == 1))];
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
        model = simcafit(x_in,y_in,comp,pret_type,assign_method);
        out = simcapred(x_out,model);
        assigned_class(find(in == 0)) = out.class_pred;
    end
end

class_param = calc_class_param(assigned_class',class);

cv.class_pred = assigned_class';
cv.class_param = class_param;
cv.settings.cv_groups = cv_groups;
cv.settings.cv_type = cv_type;
cv.settings.num_comp = comp;
cv.settings.scal = pret_type;