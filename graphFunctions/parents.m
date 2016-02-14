function ps = parents(adj_mat, ilist)
% PARENTS Return the list of parents of node i
% ps = parents(adj_mat, i)

ps = [];
for i = ilist
    ps = [ps find(adj_mat(:,i))'];
end

ps = unique(ps);
