function Tree = BackUp(Tree, childToParent, t, index)

% Starts with a pre decision state

while (Tree(index).t > t)
    Tree(index).N = Tree(index).N + 1;

    % Move backwards to previous post-decision state
    index = childToParent(index);
    Tree(index).N_x = Tree(index).N_x + 1;
    Tree(index).V_x = computePostDecisionValue();
    V_x = Tree(index).V_x;
   
    
    % Move backwards to previous pre-decision state
    index = childToParent(index);
    if Tree(index).t == t
        Tree(index).N = Tree(index).N + 1;
    end
    Tree(index).V_x = updateValue(V_x);
end

% Will have to update tehis function if the sampling distribution is
% not uniform
    function V_x = computePostDecisionValue() % Computes the value for taking a decision
        % Assuming the importance sampling uses a uniform distribution

        % Expectation taken over all explored exogenous events
        probSum = 0;
        scaling = zeros(1, length(Tree(index).Omega_e));
        for i = 1: length(Tree(index).Omega_e)
            observation = Tree(index).Omega_e(i);
            scaling(i) = Tree(index).pmf(observation)...
                /(1/length(Tree(index).Omega_e));
            probSum = probSum + Tree(index).pmf(observation);
            values(i) = Tree(Tree(index).actionChildren(observation)).V_x;

        end
        E_g     = sum(scaling.*values)/length(Tree(index).Omega_e);
        V_x     = E_g/probSum;
    end

    function V = updateValue(V_x) % Updates value at pre-decision state
        Delta = V_x;
        V = Tree(index).V_x + (Delta - Tree(index).V_x)/Tree(index).N;
    end

end