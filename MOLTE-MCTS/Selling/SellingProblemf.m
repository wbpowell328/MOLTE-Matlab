function results = SellingProblemf(type, rolloutPolicy, d_thr, e_thr, alpha, budget)

% SellingProblemf Solves a selling problem using MCTS
%   Inputs
%     type          = D for deterministic, S for stochastic, defaults
%                     to stochastic
%     rolloutPolicy = the rollout policy chosen (a function handle)
%     d_thr         = threshold for decisions examined
%     e_thr         = threshold for observations sampled
%     alpha         = parameter for UCT
%     budget        = number of iterations for each call of the MCTS within
%                     the problem
%   Outputs
%     takeAction    = the decision vector (solution)
%     p             = the prices at time t (used for checking)
%     endT          = time taken to solve the whole problem

startT = tic;
t_final = 8; % Horizon
inventory = 30;

d = 12*ones(t_final, 1);
p = 5*abs(sin(linspace(pi/4,2*pi,t_final)));
% p = [1 7 2 6 4]';

% Ensures p is a column vector
p = fixDimensions(p);
d = fixDimensions(d);

totT = 0;
forecast = struct('demand', d, 'price', p');

t = 1;

takeAction = zeros(1, t_final);
subplot(2,2,4)

[optVal, optDecision] = DeterministicIPSelling(inventory, d, p', 1e5);


for i = 1:t_final-1
    clear Tree;
    clear childToParent;
    
    Tree(1) = State_Selling();
    Tree(1) = Tree(1).PreDecisionState(0:2:min(d(i), inventory), 1);
    Tree(1).inv = inventory;
    Tree(1).t = t;
    Tree(1).N = 0;
    childToParent(1) = 0;
    iterCount = 1;
    % Makes one decision
    start = tic;
    rootValue = zeros(1, budget);
    while iterCount <= budget
        [Tree, childToParent, index] = TreePolicy_Selling(Tree, t_final, d_thr, e_thr, alpha, 1, childToParent, forecast, type);
        Tree(index) = SimPolicy_Selling(Tree(index), t_final, forecast, 1e5, rolloutPolicy);
        Tree = BackUp(Tree, childToParent, t, index);
        rootValue(iterCount) = Tree(1).V_x;
        iterCount = iterCount + 1;
    end
    subplot(2,2,4)
    hold on;
    plot(1:budget, rootValue);
    
    % Extract first decision
    rootChildren = Tree(1).children;
    bestAction = 2;
    for j = 1:length(rootChildren)
        if Tree(bestAction).V_x < Tree(rootChildren(j)).V_x
            bestAction = rootChildren(j);
        end
    end
    endT = toc(start);
    totT = totT + endT;
    
    takeAction(t) = Tree(bestAction).action;
    
    % Update inventory
    inventory = inventory - takeAction(t);
    fprintf('Updated inventory: %d\n', inventory);
    if inventory <= 0
        break;
    end
    t = t + 1;
end

if inventory > 0
    takeAction(t_final) = min(d(end), inventory);
end

endT = toc(startT);

%% Storing results
bestSolution    = struct('decision', takeAction, 'value', takeAction*p);
optSolution     = struct('decision', (optDecision)', 'value', optVal);
results         = struct('bestSolution', bestSolution, 'optSolution', optSolution, 'time', endT);

%% Plotting
% Decisions
subplot(2,2,1)
hold on;
plot(1:t_final, takeAction, '-x');
plot(1:t_final, optDecision, '-o');
title('Decisions');
xlabel('Time');
ylabel('Decision');
legend('MCTS', 'Optimal', 'Location', 'best');

bestRev = takeAction'.*p;
optRev = optDecision.*p;

% Reward
subplot(2,2,2)
hold on;
plot(1:t_final, bestRev, '-x');
plot(1:t_final, optRev, '-o');
title('Revenue');
xlabel('Time');
ylabel('Revenue');
legend('MCTS', 'Optimal', 'Location', 'best');

% Cumulative Revenue
subplot(2,2,3)
hold on;
plot(1:t_final, cumsum(bestRev), '-x');
plot(1:t_final, cumsum(optRev), '-o');
title('Cumulative Revenue');
xlabel('Time');
ylabel('Revenue');
legend('MCTS', 'Optimal', 'Location', 'best');
% 
subplot(2,2,4)
title('Convergence');
xlabel('Iteration');
legend(cellstr(num2str((1:i)'))', 'Location', 'best');


% Flips a row vector to column vector
    function x = fixDimensions(x)
        [nr, nc] = size(x);
        if nc > nr
            x = x';
        end
    end


    
end

