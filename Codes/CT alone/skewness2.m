function s = skewness2(x,flag,dim)
%SKEWNESS Skewness.
%   S = SKEWNESS(X) returns the sample skewness of the values in X.  For a
%   vector input, S is the third central moment of X, divided by the cube
%   of its standard deviation.  For a matrix input, S is a row vector
%   containing the sample skewness of each column of X.  For N-D arrays,
%   SKEWNESS operates along the first non-singleton dimension.
%
%   SKEWNESS(X,0) adjusts the skewness for bias.  SKEWNESS(X,1) is the same
%   as SKEWNESS(X), and does not adjust for bias.
%
%   SKEWNESS(X,FLAG,'all') is the skewness of all the elements of X.
%
%   SKEWNESS(X,FLAG,DIM) takes the skewness along dimension DIM of X.
%
%   SKEWNESS(X,FLAG,VECDIM) finds the skewness of the elements of X based 
%   on the dimensions specified in the vector VECDIM.
%
%   SKEWNESS treats NaNs as missing values, and removes them.
%
%   See also MEAN, MOMENT, STD, VAR, KURTOSIS.

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
    if isequal(x,[]), s = NaN('like',x); return; end

    % Figure out which dimension nanmean will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Center X, compute its third and second moments, and compute the
% uncorrected skewness.
% x0 = x - repmat(nanmean(x,dim), t);
x0 = x - nanmean(x,dim);
s2 = nanmean(x0.^2,dim); % this is the biased variance estimator
m3 = nanmean(x0.^3,dim);
s = m3 ./ (s2.^(1.5)+1e-10);
% Bias correct the skewness.
if flag == 0
    n = sum(~isnan(x),dim);
    n(n<3) = NaN; % bias correction is not defined for n < 3.
    s = s .* sqrt((n-1)./n) .* n./(n-2);
end
end
