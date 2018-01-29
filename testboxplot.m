function [Q1 Q3 outlier outlier_num] = testboxplot(x, x_date, k, tukey, bplot)
%
% The 'testboxplot' identifies outliers that lie 'k' times over the
% interquantile rante (Q3-Q1). A boxplot using the standard Matlab function
% 'boxplot' is also provided by default.
%
% Data in 'x' are organized so that columns are the time series and rows
% are the time intervals. All series contain the same number of
% observations.
%
% 'x_date' is a column vector (cell array) containing the dates of the
% corresponding data in 'x'.
%
% The variable 'k' is a multiplier of the interquantile range (Q3-Q1). The
% default value is 2.
%
% The function also provides an option ('tukey') to use the quantile
% computation suggested by Tukey. It can take the default value 0 (not
% applied) or 1 (apply Tukey approach).
%
% The variable 'bplot' indicates if a boxplot is provided (=1, default) or
% omitted (=0).
%
% [Q1 Q3 outlier outlier_num] = testboxplot(...) returns the following information:
% Q1            - 25% quantile
% Q3            - 75% quantile
% outlier       - cell array specifying the date and the series number (column
%                 number in 'x') where the potential outliers are situated.
% outlier_num   - matrix providing row and column numbers of the values in
%                 'x' considered as potential outliers.
%
% Created by Francisco Augusto Alcaraz Garcia
%            alcaraz_garcia@yahoo.com
%
% References:
%
% 1) B. Iglewicz; D.C. Hoaglin (1993). How to Detect and Handle Outliers.
% ASQC Basic References in Quality Control, vol. 16, Wisconsin, US.
%
% 2) J.W. Tukey (1977). Exploratory Data Analysis. Addison Wesley.


% Check number of input arguements
if (nargin < 2) || (nargin > 5)
    error('Requires two to five input arguments.')
end

% Define default values
if nargin == 2,
    k = 2;
    tukey = 0;
    bplot = 1;
elseif nargin == 3,
    tukey = 0;
    bplot = 1;
elseif nargin == 4,
    bplot = 1;
end

% Check for validity of inputs
if ~isnumeric(x) || ~iscellstr(x_date),
    error('Input x must be a numeric array and x_date must be a string table.')
end

[n, c] = size(x);

if tukey == 0,
    Q1 = quantile(x, 0.25);
    Q3 = quantile(x, 0.75);
else
    [xsort, ix] = sort(x);
    f = ((n+1)/2+1)/2;
    if isinteger(f),
        Q1 = xsort(f,:);
        Q3 = xsort(end+1-f,:);
    else
        Q1 = (xsort(fix(f),:) + xsort(fix(f)+1,:))/2;
        Q3 = (xsort(end+1-fix(f),:) + xsort(end+1-fix(f)+1,:))/2;
    end
end

[i1,j1] = find(x < repmat(Q1-k*(Q3-Q1),n,1));
[i2,j2] = find(x > repmat(Q3+k*(Q3-Q1),n,1));

if (isempty(i1)+isempty(i2)) == 0,
    outlier = [x_date(i1) cellstr(strcat('Series', num2str(j1)));
    x_date(i2) cellstr(strcat('Series', num2str(j2)))];
    outlier_num = [i1 j1; i2 j2];
elseif isempty(i1),
    outlier = [x_date(i2) cellstr(strcat('Series', num2str(j2)))];
    outlier_num = [i2 j2];
elseif isempty(i2),
    outlier = [x_date(i1) cellstr(strcat('Series', num2str(j1)))];
    outlier_num = [i1 j1];
else
    outlier = ('No outliers have been identified!');
    outlier_num = ('No outliers have been identified!');
end

if bplot == 1,
    boxplot(x, 'notch', 'on', 'whisker', k, 'symbol', 'r.')
end

