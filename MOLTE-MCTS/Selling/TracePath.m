function path = TracePath(Tree, index, t_final, path)
if Tree(index).t == t_final
    return;
end


currentChildren = Tree(index).children;
for i = 1:length(currentChildren)
    childIndex = currentChildren(i);
    if i == 1
        maxCount = childIndex;
    elseif Tree(childIndex).V_x > Tree(maxCount).V_x
        maxCount = childIndex;
    end
end
    path(length(path)+1) = Tree(maxCount).action;
%     fprintf('Node #: %d, City: %d', maxCount, Tree(maxCount).action);

    path = TracePath(Tree, Tree(maxCount).children, t_final, path); % Move to next post decision state
end