function model = dafit(X,class,class_prob,method,num_comp,scal)

% Discriminant Analysis fitting
%
% model = dafit(X,class,class_prob,method,num_comp,scal)
%
% input:
% X             dataset [samples x variables]
% class         class vector [samples x 1]
% class_prob    if prob = 1 equal probability, if prob = 2 proportional prob.
% method        'linear' or 'quadratic'
% optional:
% num_comp      DA is calculated on the first num_comp principal components
% scal          scaling method when PCA is calculated
%               if scal = 'none' -> no scaling
%               if scal = 'cent' -> centering
%               if scal = 'auto' -> autoscaling
%
% output:
% model is a structure conyaining
% class_calc    calculated class vector [n x 1]
% class_param   classification parameters
% settings      model settings
% L             loadings on canonical variables [p x G-1], only for LDA
% Lstd          standardized loadings on canonical variables [p x G-1], only for LDA
% S             scores on canonical variables [n x G-1], only for LDA
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

[nobj,nvar] = size(X);
nclass = max(class);

if class_prob == 2
    for g = 1:nclass 
        obj_cla(g)  = length(find(class == g));
    end
    prob = obj_cla/nobj;
end

if nargin == 6
    modelpca = pca_model(X,num_comp,scal,0);
    Xtrain = modelpca.T;
else
    Xtrain = X;
    modelpca = NaN;
    num_comp = 0;
end

% if linear and not with PCs check for pooled estimate of covariance
doit = 1;
if strcmp('linear',method) & nargin < 6
    doit = pec(X,class);
end

if doit
    % fitting
    if class_prob == 1
        [class_calc] = classify(Xtrain,Xtrain,class,method);
    else
        [class_calc] = classify(Xtrain,Xtrain,class,method,prob);
    end
    % calculates canonical variables for lda
    if strcmp('linear',method)
        class_unfold = zeros(size(Xtrain,1),max(class));
        for g=1:max(class)
            class_unfold(find(class==g),g) = 1;
        end
        [L,B,r,S,V] = canoncorr(Xtrain,class_unfold);
        for k=1:size(L,1);
            for j=1:size(L,2)
                Lstd(k,j) = L(k,j)*std(Xtrain(:,k));
            end
        end
    end
else
    class_calc = ones(size(X,1),1);
    L = zeros(size(X,2),1);
    Lstd = zeros(size(X,2),1);
    S = zeros(size(X,1),1);
end

class_param = calc_class_param(class_calc,class);
settings.class_prob = class_prob;
settings.method = method;
settings.Xtrain = Xtrain;
settings.modelpca = modelpca;
settings.num_comp = num_comp;
settings.class_true = class;

model.class_calc  = class_calc;
model.class_param = class_param;
model.settings = settings;

if strcmp('linear',method)
    model.L = L;
    model.Lstd = Lstd;
    model.S = S;
end

% -------------------------------------------------------------------------
function doit = pec(X,class)

for g = 1:max(class)
    gmeans(g,:) = mean(X(find(class == g),:));
end
% Pooled estimate of covariance
[Q,R] = qr(X - gmeans(class,:), 0);
R = R / sqrt(size(X,1) - max(class)); % SigmaHat = R'*R
s = svd(R);
if any(s <= eps^(3/4)*max(s))
    doit = 0;
    disp('The pooled covariance matrix of training samples must be positive definite. model not calculated');
else
    doit = 1;
end