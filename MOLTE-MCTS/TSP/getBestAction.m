function bestAction = getBestAction(Tree, children)
    bestAction = children(1);
    for i = 1:length(children)
        if Tree(children(i)).V_x <= Tree(bestAction).V_x
            bestAction = Tree(children(i)).action;
        end
    end
end