function [X,Sk_sorted,F] = sindyAlgorithm(A,b,lambda)

% ===================================================================
% Description: 
%   This function solves the linear system Ax=b using the 
%   following algorithm: 
%       x^{0}   = A\b
%       S^{k}   = {j: |x^{k}_{j}| >= lambda}
%       x^{k+1} = argmin || Ax-b ||_2 s.t. supp(x) in S^{k}
%
% Inputs:
%   A = matrix of size m*n
%   b = vector of size m*1
%   lambda = thresholding parameter
%
% Outputs:
%   X = the solution at each iteration
%   Sk_sorted = the ordered support set of each x^k
%   F = the value of the objective function at each step
%
% References:
%   [1] Steven L. Brunton, Joshua L. Proctor, and J. Nathan Kutz. 
%       Discovering governing equations from data by sparse 
%       identification of nonlinear dynamical systems. Proceedings of
%       the National Academy of Sciences, 113(15):3932-3937, 2016.
%   [2] Linan Zhang and Hayden Schaeffer. 
%       On the Convergence of the SINDy Algorithm. Multiscale Modeling
%       & Simulation 17(3), 948–972, 2019.
%
% Authors: Linan Zhang and Hayden Schaeffer
% Date: May 16, 2018
% ===================================================================

% Set iteration parameters.
n = size(A,2); MaxIt = n; % maximum number of iterations needed
err = 1; % flag for error
k = 2; % iteration index

% Initialize x^{k} and S^{k} for each k.
X = zeros(n,MaxIt); % Xk(:,k) = x^{k-1}
Sk = zeros(n,MaxIt); % Sk(:,k) = S^{k-1}
% Sort Sk such that |x^{k}_{S^{k}(j+1)}| >= |x^{k}_{S^{k}(j)}|
Sk_sorted = zeros(n,MaxIt);

% Perform the initial step.
x = A\b; % x = x^{0}
X(:,1) = x;
S = find( abs(x)>=lambda ); % S = S^{0}
% Sort S such that |x_{S(j+1)}| >= |x_{S(j)}|.
[~,ind] = sort(abs(x(S)),'descend');
Sk(1:length(ind),1) = S;
Sk_sorted(1:length(ind),1) = S(ind);

% Perform the iterative step.
while k<=MaxIt && err
    
    % xnew = argmin ||Ay-b||_2 subject to supp(y) in S
    y = A(:,S)\b;
    x = zeros(n,1);
    x(S) = y; % x = x^{k-1}
    X(:,k) = x;
    S = find( abs(x)>=lambda ); % S = S^{k-1}
    % Sort S such that |x_{S(j+1)}| >= |x_{S(j)}|.
    [~,ind] = sort(abs(x(S)),'descend');
    Sk(1:length(ind),k) = S;
    Sk_sorted(1:length(ind),k) = S(ind);

    % Stopping criterion: S^{k-1} = S^{k-2}.
    err = ~isequal(Sk(:,k-1), Sk(:,k));
    k = k+1;
    
end

Sk_sorted(:,k:end) = [];
Sk_sorted = Sk_sorted(any(Sk_sorted,2),:);
X(:,k:end) = [];

% Compute the objective function.
k = size(X,2); % number of iterations (including the 0th step)
normA = norm(A,2); % l2 norm of A
F = zeros(k,1);
for ii=1:k
    F(ii) = obj(X(:,ii));
end
    function f = obj(x)
        f = norm(A*x-b,'fro')^2 / (normA^2) + lambda^2 * nnz(x);
    end

end
