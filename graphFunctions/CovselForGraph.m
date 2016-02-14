function [Z, runtime] = CovselForGraph(S, lambda, rho, alpha,Domain,TOL)

% Modification of the covsel function found in the paper:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
% Only cares about solving for the graph
% Solves the following problem via ADMM:
%
%   minimize  trace(S*X) - log det X + lambda*||X_Domc||_1
%
%   subject to X_Dom0 = 0.

% where S is the empirical covariance of the data.  Thus, we only
% enforce sparsity on entries in Domc and know that the entries in
% Dom0 are zero.
%
% Domain specifies Domc and Dom0.
% Domain(i,j) = 1 => (i,j) in Dom0
% Domain(i,j) = 0 => (i,j) in Domc
% Domain(i,j) = 2 => X(i,j) is non-zero.

% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html

% Global constants and defaults

time0 = cputime;
MAX_ITER = 1000;

% Data preprocessing
n = size(S,1);

% ADMM solver
if ~isempty(Domain)
    Domain = setdiag(Domain,2);
else
    Domain = speye(n,n)*2;
end

Z = zeros(n);
U = zeros(n);

Dom = (Domain == 2); %X(i,j) not zero
Dom0 = (Domain == 1); % X(i,j) zero
Domc = (Domain == 0); % X(i,j) not known, depends on lambda

if isempty(Domc) % all sparsity pattern is known
    Z = Domain - 1;
end

if sum(Domc(:)) == 0
    CheckDomain = Dom+Dom0;
else
    CheckDomain = Domc;
end

for k = 1:MAX_ITER

    % x-update
    [Q,L] = eig(rho*(Z - U) - S);
    Q = real(Q); L = real(L);
    es = diag(L);
    xi = (es + sqrt(es.^2 + 4*rho))./(2*rho);
    X = Q*diag(xi)*Q';

    % z-update
    Zold = Z;
    X_hat = alpha*X + (1 - alpha)*Zold;
    Z(Domc) = shrinkage(X_hat(Domc) + U(Domc), lambda/rho);
    Z(Dom0) = 0;
    Z(Dom) = shrinkage(X_hat(Dom) + U(Dom), 0/rho);

    % u-update
    U = U + X_hat - Z;

    % check convergence
    ss_norm = norm(X(CheckDomain) - Z(CheckDomain)) / length(CheckDomain);
%    ss_norm = norm(X - Z);

    if (ss_norm < TOL)
        break;
    end

end

runtime = cputime - time0;

end

function y = shrinkage(a, kappa)
y = max(0, a-kappa) - max(0, -a-kappa);
end
