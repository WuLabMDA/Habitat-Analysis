function k = kurtosis2(x,flag,dim)
%KURTOSIS Kurtosis.
%   K = KURTOSIS(X) returns the sample kurtosis of the values in X.  For a
%   vector input, K is the fourth central moment of X, divided by fourth
%   power of its standard deviation.  For a matrix input, K is a row vector
%   containing the sample kurtosis of each column of X.  For N-D arrays,
%   KURTOSIS operates along the first non-singleton dimension.
%
%   KURTOSIS(X,0) adjusts the kurtosis for bias.  KURTOSIS(X,1) is the same
%   as KURTOSIS(X), and does not adjust for bias.
%
%   KURTOSIS(X,FLAG,'all') is the kurtosis of all the elements of X.
%
%   KURTOSIS(X,FLAG,DIM) takes the kurtosis along dimension DIM of X.
%
%   KURTOSIS(X,FLAG,VECDIM) finds the kurtosis of the elements of X based
%   on the dimensions specified in the vector VECDIM.
%
%   KURTOSIS treats NaNs as missing values, and removes them.
%
%   See also MEAN, MOMENT, STD, VAR, SKEWNESS.

%   Copyright 1993-2018 The MathWorks, Inc.


if nargin < 2 || isempty(flag)
    flag = 1;
end

% Validate flag
if ~(isequal(flag,0) || isequal(flag,1) || isempty(flag))
    error(message('stats:trimmean:BadFlagReduction'));
end

if nargin < 3 || isempty(dim)
    % The output size for [] is a special case, handle it here.
    if isequal(x,[]), k = NaN('like',x); return; end

    % Figure out which dimension nanmean will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Center X, compute its fourth and second moments, and compute the
% uncorrected kurtosis.
x0 = x - nanmean(x,dim);
s2 = nanmean(x0.^2,dim); % this is the biased variance estimator
m4 = nanmean(x0.^4,dim);
k = m4 ./ (s2.^2+1e-10);

% Bias correct the kurtosis.
if flag == 0
    n = sum(~isnan(x),dim);
    n(n<4) = NaN; % bias correction is not defined for n < 4.
    k = ((n+1).*k - 3.*(n-1)) .* (n-1)./((n-2).*(n-3)) + 3;
end
end
