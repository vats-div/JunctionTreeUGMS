function H = MuGScreeningGraphs(X,K,m)

% MuG Screening for Learning Graphs

if nargin == 1
    K = 5;
    m = 2;
end

p = size(X,2);
H = zeros(p,p);

% Run MuG for each node

for k = 1:p
    y = X(:,k);
    XX = X;
    XX(:,k) = [];
    Sbar = MuGScreeningGroupLasso(y,XX,K,m);
    IndexAbovek = (Sbar >= k);
    Sbar(IndexAbovek) = Sbar(IndexAbovek) + 1;
    H(k,Sbar) = H(k,Sbar) + 1;
    H(Sbar,k) = H(Sbar,k) + 1;
end

H = triu(H) + triu(H)';
H = double(H == 2);

end