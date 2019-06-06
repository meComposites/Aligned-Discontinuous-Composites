function y = SPlinspace2(d1, d2, n)
%LINSPACE2 Linearly spaced vector.
%   LINSPACE(X1, X2) generates a matrix of 100 linearly
%   equally spaced points between X1 and X2.
%
%   LINSPACE(X1, X2, N) generates N points between X1 and X2.
%   For N < 2, LINSPACE returns X2.
%
%   Class support for inputs X1,X2:
%      float: double, single
%
%   Example 1: Vector starting and ending points
%       d1 = [0:0.1:10];
%       d2 = [10:-0.1:0];
%       M = linspace2(d1,d2);
%       imagesc(M);
%
%   Example 2: Mixing Vectors and Scalars
%       d1 = 0;
%       d2 = [10:-0.1:0];
%       M = linspace2(d1,d2);
%       imagesc(M);

%   Jonathan Sullivan
%   MIT Lincoln Laboratory
%   jonathan.sullivan@ll.mit.edu
%   Original 10-Feb-2011

if nargin == 2
    n = 100;
end

n = double(n);
ld1 = length(d1);
ld2 = length(d2);

if ld1 == 1 && ld2 > 1
    d1 = repmat(d1,ld2,1);
elseif ld2 == 1 && ld1 > 1
    d2 = repmat(d2,ld1,1);
end

d1 = d1(:);
d2 = d2(:);

d1a = repmat(d1,1,n-1);
y = [d1a+(d2-d1)*(0:n-2)/(floor(n)-1) d2];