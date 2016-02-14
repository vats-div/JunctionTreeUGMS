function ss = randomSign(m,n)

% create a random mxn matrix of +1,-1

if nargin == 1
    n = 1;
end

ss = double(rand(m,n) > 0.5);

ss(ss == 0) = -1;

end