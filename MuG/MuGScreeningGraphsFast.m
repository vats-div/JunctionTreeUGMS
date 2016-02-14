function H = MuGScreeningGraphsFast(X,K,m,kE)

% MuG Screening for Learning Graphs

if nargin == 1
    K = 5;
    m = 2;
end

if nargin == 3
    kE = size(X,1);
end

p = size(X,2);
H = zeros(p,p);

% Run MuG for each node

for k = 1:p
    y = X(:,k);
    XX = X;
    XX(:,k) = [];
    display(['node k = ' num2str(k)]);
    Sbar = MuGScreeningGroupLassoFast(y,XX,K,m,kE);
    IndexAbovek = (Sbar >= k);
    Sbar(IndexAbovek) = Sbar(IndexAbovek) + 1;
    H(k,Sbar) = H(k,Sbar) + 1;
    H(Sbar,k) = H(Sbar,k) + 1;
end

H = triu(H) + triu(H)';
H = double(H == 2);

end