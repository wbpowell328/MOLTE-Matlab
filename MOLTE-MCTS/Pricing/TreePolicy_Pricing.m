function [Tree, childToParent, index] = TreePolicy_Pricing(Tree, t_final, d_thr, e_thr, alpha, MVar, childToParent, theta_truth)
numNodes = length(Tree);
index = 1;
if isempty(childToParent)
    childToParent(index) = index;
end

% Traverse down the tree till leaf node
while (Tree(index).t < t_final)
    %% Expand from pre decision
    
    % if I have less decisions used at this state then I wanted, I do an expansion
    if (length(Tree(index).X_e) < d_thr && ~isempty(Tree(index).X_u))

        % Randomly select unused action, x is actual decision
        x = Tree(index).X_u(randi(length(Tree(index).X_u)));
        
        % Maps nodes to an action taken
        UpdateActionChildren(x);

        % Create post decision state
        ExpandAction();
        
        % Points to newest node
        index = numNodes;
    else
        % Move to post decision state as selected by UCT
        x = UCT(alpha, 'max');
        
        % Points to node chosen by UCT
        index = Tree(index).actionChildren(Tree(index).X_e(x));
    end
    
    %% Expand from post decision
    % If condition: explored observations >= e_thr and there are still
    % unexplored observations
    if length(Tree(index).Omega_e) < e_thr && ~isempty(Tree(index).Omega_u)
        
        % Choose a random exogenous event with uniform distribution
        w = Tree(index).Omega_u(randi(length(Tree(index).Omega_u)));

        % Maps nodes to observation chosen
        UpdateActionChildren(w);

        % Realizes observation and create new pre-decision state
        ExpandObservation();
        
        % Points to newest node
        index = length(Tree);
        
        % end
        return
    else
        
        % Choose random exogenous event that has been explored
        w = Tree(index).Omega_e(randi(length(Tree(index).Omega_e)));
        
        % Points to node
        index = Tree(index).actionChildren(w);
    end
end

%% Expands from pre-decision to post-decision
    function ExpandAction()
        numNodes = numNodes + 1;
        Tree(numNodes) = State_Pricing();
        Tree(numNodes) = TransitionPre2Post_Pricing(Tree(index), Tree(numNodes));                     
        [Tree(numNodes), Omega] = generateObservations_Pricing(Tree(numNodes), theta_truth, MVar, x);
        Tree(numNodes) = Tree(numNodes).PostDecisionState(Omega, numNodes);
        Tree(numNodes).action = x;
        Tree(numNodes).t = Tree(index).t;
        Tree(index).children(length(Tree(index).children)+1) = numNodes;
        childToParent(numNodes) = index;
        UpdateActionSet(); % Update X_e and X_u
    end


%% Expands from post-decision to pre-decision
    function ExpandObservation()
        numNodes = numNodes + 1;
        Tree(numNodes) = State_Pricing();
        decisions = generateDecisions_Pricing();
        Tree(numNodes) = Tree(numNodes).PreDecisionState(decisions, numNodes);
        Tree(numNodes) = TransitionPost2Pre_Pricing(Tree(index), Tree(numNodes), w, 0.5);    % Realizes profit        
        Tree(numNodes).realization = w;
        Tree(index).children(length(Tree(index).children)+1) = numNodes;
        Tree(numNodes).t = Tree(index).t + 1;
        childToParent(numNodes) = index;
        UpdateObservationSet();
    end

    function UpdateActionSet()
        Tree(index).X_e(length(Tree(index).X_e) + 1) = x;
        Tree(index).X_u = setdiff(Tree(index).X_u, x);
    end

    function UpdateObservationSet()
        Tree(index).Omega_e(length(Tree(index).Omega_e) + 1) = w;
        Tree(index).Omega_u = setdiff(Tree(index).Omega_u, w);
    end

    function x = UCT(alpha, type)
        value= zeros(1,length(Tree(index).X_e));
        w    = zeros(1,length(Tree(index).X_e));
        
        for i = 1:length(Tree(index).X_e)
            decision = Tree(index).X_e(i);
            postDecisionStateInd = Tree(index).actionChildren(decision);
            w(i) = sqrt(log(Tree(index).N)/Tree(postDecisionStateInd).N_x);
            value(i) = Tree(postDecisionStateInd).V_x;
        end
        
        % let cost be zero for this problem
        cost = 0;
        if strcmp(type, 'max')
            [~, x] = max((cost + value) + alpha .* w); % Maximization
        else
            if ~strcmp(type, 'min')
                fprintf('UCT type: %s invalid. Defaulting to minimization', type);
            end
            [~, x] = max(-(cost + value) + alpha .* w); % Minimization
        end
    end

    function UpdateActionChildren(arg)
        newKeys = [Tree(index).actionChildren.keys, arg];
        newVals = [Tree(index).actionChildren.values, numNodes + 1];
        
        % Links decisions (outcomes respectively) to post (pre respectively)
        % decision states
        Tree(index).actionChildren = containers.Map(newKeys, newVals);
    end
end