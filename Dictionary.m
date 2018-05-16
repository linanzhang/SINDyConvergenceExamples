function [D,L,D_cnorm] = Dictionary(u,p,varargin)

% ===================================================================
% Description:
%   The function outputs the Monomial Dictionary Matrix.
%
% Inputs:
%   u = the state of a nonlinear system, i.e. given data
%   p = maximum degree of polynomials in the dictionary
%
% Optional input:
%   doNormalization = flag for normalizing the columns of D
%
% Outputs:
%   D = dictionary of monomials of degree at most p
%   L = legend of D
%   D_cnorm = the norm of each column of the unscaled D.
%
% Example: 
%   Inputs: u=[u1,u2,u3], p=2
%   Output: D = [1, u1, u2, u3, u1^2, u1u2, u1u3, u2^2, u2u3, u3^2]
%
% Remark: 
%   The function can be modified to output the Trigonometric
%   Dictionary Matrix.
%
% Reference:
%   https://www.math.cmu.edu/CNA/Publications/publications2018/papers/18-CNA-013.pdf
%
% Authors: Linan Zhang and Hayden Schaeffer
% Date: May 16, 2018
% ===================================================================

% Compute the size of u.
%   m = number of time steps
%   n = number of states.
[m,n] = size(u);
% Augment u to [1;u].
u1 = [ones(m,1) u];

% List the grid points.
l1 = cell(n+1,1);
if n==2
    % If only 2 states, use "u" and "v" to denote the states.
    l1{2} = 'u';
    l1{3} = 'v';
else
    % Otherwise, use "u1", "u2", etc. to denote the states.
    for ii=1:n
        l1{ii+1} = strcat('u',num2str(ii));
    end
end

% Find all polynomials of degree at most p.
C = 1:(n+1);
for ii=1:p-1
    C = combvec(C,1:(n+1));
end
% Sort each column of C in ascending order.
C = sort(C,1);
% Remove duplicate rows in C' (i.e. duplicate columns in C).
C = unique(C','rows');
% number of rows of C = number of columns of D
% number of columns of C = p
nD = size(C,1);

% Construct D.
D = ones(m,nD);
L = cell(nD,1);
for ii=1:nD
    % Multiply the corresponding columns to make a polynomial of 
    % degree at most p.
    for jj = 1:p
        D(:,ii) = D(:,ii) .* u1(:,C(ii,jj));
        L{ii} = strcat(L{ii},l1{C(ii,jj)});
    end
end
L{1} = num2str(1);

% Normalize D.
if nargin==2
    doNormalization = 0;
else 
    doNormalization = varargin{1};
end
if doNormalization
    % Compute the norm of each column.
    D_cnorm = sqrt(sum(D.^2,1)); 
    D_cnorm = D_cnorm(:);
    % Normalize columns of D.
    D = normc(D);
end
    