startT = tic;
d_thr = 10;
e_thr = 5;
alpha = .5;
budget = 20;
problemType = 'D';
rolloutPolicty = str2func('DeterministicTSP');
% dpath = {'EqualByValueDoubleSet.jar'};
% javaclasspath(dpath);
n_cities = 4;
costMat = [0 2 9 10; 1 0 6 4; 15 7 0 8; 6 3 12 0];
% costMat = 10*rand(n_cities, n_cities);
% costMat = [0 12 11 16; 15 0 15 10; 8 14 0 18; 9 11 17 0];
% n_cities = 4;
% load('costMat6City.mat');
t_final = n_cities;
bestAction = zeros(n_cities+1, 1);
bestAction([1, end],1) = 1;
accCost = 0;

visitedCities = 1;
for t = 1:t_final-1
    
    clear Tree;
    clear childToParent;
    Tree(1) = State_TSP();
    Tree(1).city = bestAction(t);
    
    
    decisions = setdiff((1:n_cities)', bestAction);
    Tree(1) = Tree(1).PreDecisionState(decisions,1);
    Tree(1).t = t;
    
    Tree(1).visited = visitedCities;
    Tree(1).accCost = accCost;
    fprintf('Accumulated Cost: %d\n', accCost);
    childToParent(1) = 0;
    
    iter = 1;
    rootNodeVal = zeros(budget, 1);
    while (iter <= budget)
        [Tree, childToParent, index] = TreePolicy_TSP(Tree, t_final, d_thr, e_thr, alpha, 1, childToParent, costMat, problemType);
        Tree(index) = SimPolicy_TSP(rolloutPolicy, Tree(index), 1, n_cities, costMat);
        fprintf('Tree(%d).V_x = %d\n', index, Tree(index).V_x);
        Tree = BackUp(Tree, childToParent, t, index);
        rootNodeVal(iter) = Tree(1).V_x; 
%         fprintf('Root node val: %d\n', Tree(1).V_x);
        iter = iter + 1;
    end
    
    bestAction(t+1) = getBestAction(Tree, Tree(1).children);
    accCost = accCost + costMat(bestAction(t), bestAction(t+1));
    visitedCities(length(visitedCities)+1) = bestAction(t+1);
end

endT = toc(startT);
accCost = accCost + costMat(bestAction(end-1), bestAction(end));
% Optimal Solution
[optPath, val, ~] = TSPDP(1:n_cities, 1, costMat, 1000);