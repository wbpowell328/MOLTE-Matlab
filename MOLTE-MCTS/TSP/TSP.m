
clear all

startT = tic;

dpath = {'EqualByValueDoubleSet.jar'};
javaclasspath(dpath);

costMat = [0 2 9 10; 1 0 6 4; 15 7 0 8; 6 3 12 0];

n_cities = 4;

d_thr = n_cities;
e_thr = 3;
alpha = 0.5;
type = 'D';

t_final = n_cities;
bestAction = zeros(n_cities+1, 1);
bestAction([1, end],1) = 1;
accCost = 0;
for t = 1:t_final-1
    
    Tree(1) = State_TSP();
    if t == 1
        Tree(1).city = 1;
    else
        Tree(1).city = bestAction(t);
    end
    
    decisions = setdiff((1:n_cities)', bestAction);
    Tree(1) = Tree(1).PreDecisionState(decisions,1);
    Tree(1).t = t;
    
    Tree(1).visited = Tree(1).city;
    Tree(1).accCost = accCost;
    
    childToParent(1) = 0;
    iter = 0;
    index = 1;
    
    values = zeros(250,1);
    while (iter <= 300)
        [Tree, childToParent, index] = TreePolicy_TSP(Tree, t_final, d_thr, e_thr, alpha, 1, childToParent, costMat, type);
        Tree(index) = SimPolicy_TSP(Tree(index), Tree(1).city, n_cities, costMat);
        Tree = BackUp(Tree, childToParent, t, index);
        iter = iter + 1;
        values(iter) = Tree(1).V_x;
    end
    
    bestAction(t+1) = getBestAction(Tree, Tree(1).children);
    accCost = accCost + costMat(bestAction(t), bestAction(t+1));
end

endT = toc(startT);
accCost = accCost + costMat(bestAction(end-1), bestAction(end));
% Optimal Solution
[optPath, val, ~] = TSPDP(1:n_cities, 1, costMat, 1000);
