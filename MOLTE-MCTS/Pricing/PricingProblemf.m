function results = PricingProblemf(~, rolloutPolicy, d_thr, e_thr, alpha, budget, t_final, L, param, resolution, weight)

if isnan(t_final)
    t_final = 8;
end
if isnan(L)
    L = 5;
end
if isnan(resolution)
    resolution = 0.01;
end
if isnan(param)
    param = 0.5;
end
if isnan(weight)
    weight = 0.6;
end
% Pricing Problem
% t_final   = Number of experiments to conduct
% L         = Number of steps to lookahead in rollout policy (truncated to
%             than t_final) t_final if greater

t_final = t_final + 1;

X_i1=[1;2;3;4;5;6;7;8;9;10];
X=[ones(10,1),X_i1,X_i1.^2];

%Prior Estimates
%theta_0=[64;2;-.2];
theta_0=[64;4;-.4];

%Standard deviations for Theta
%theta_Sigma=[20;5;.5];
theta_Sigma=[10;1;.1];

%Pick a true value for theta
thetaT=theta_0+randn(3,1).*theta_Sigma;

%Set the true values as X*theta
truth=X*thetaT;

%Set the measurement variance
MVar=4^2;

%Run the KGCBLinR L many times for rollout
% L=7;

%Set the prior on covariance matrix, assume independence
C=diag(theta_Sigma).^2;

% Parameters
% t_final = 7;
% d_thr = 10;
% e_thr = 5;
% alpha = 1;

% weight between 0 - 1, higher weight gives bias on info_gain

% Initialize the B, C, and theta
% num_init = 2;
Xs = X(randperm(length(X_i1),3),:);
Y = Xs*thetaT + normrnd(0, sqrt(MVar), 3, 1);

B = inv(Xs'*Xs);
theta = B*Xs'*Y;
bestSol = zeros(t_final-1,1);

startT = tic;
for t = 1:t_final-1
    %     fprintf('t = %d: \n', t);
    clear Tree;
    clear childToParent;
    
    Tree(1)         = State_Pricing();
    Tree(1)         = Tree(1).PreDecisionState(1:10, 1);
    Tree(1).B       = B;
    Tree(1).theta   = theta;
    Tree(1).C       = C;
    Tree(1).t       = t;
    
    childToParent(1) = 0;
    
    for i = 1:budget
        [Tree, childToParent, index] = TreePolicy_Pricing(Tree, t_final, d_thr, e_thr, alpha, MVar, childToParent, thetaT);
        Tree(index) = SimPolicy_Pricing(rolloutPolicy, Tree(index), Tree(childToParent(index)).theta, X, MVar, truth, L, resolution, .6, param);
        Tree = BackUp(Tree, childToParent, t, index);
    end
    
    maxVal = Tree(2).V_x;
    maxNode = 2;
    children = Tree(1).children;
    for i = 1:length(children)
        currNode = Tree(children(i));
        if maxVal < currNode.V_x
            maxVal = currNode.V_x;
            maxNode = children(i);
        end
    end
    %     fprintf('Action: %d, Value: %d, Count: %d\n', Tree(maxNode).action, Tree(maxNode).V_x, Tree(maxNode).N_x);
    theta = Tree(maxNode).theta;
    B = Tree(maxNode).B;
    C = Tree(maxNode).C;
    bestSol(t) = Tree(maxNode).action;
end
endT = toc(startT);

%% Comparing to KGCBLin
[thetaLast,KGCBLinSol]=KGCBLinPricing(X,theta,MVar,truth,t_final-1,C,0,1);

%% Storing Solutions
best = struct('decision', bestSol, 'value', theta);
KGCBLin = struct('decision', KGCBLinSol, 'value', thetaLast);
results = struct('MCTS', best, 'KGCBLin', KGCBLin, 'truth', thetaT, 'time', endT);

%% Plotting
% Price Curves
% figure('Name', 'Price Curves', 'NumberTitle', 'off');
subplot(2,2,1)
hold on;
plot(X_i1, X*theta);
plot(X_i1, X*thetaLast);
plot(X_i1, truth);
title('Price Curves');
xlabel('x');
ylabel('y');
legend('MCTS', 'KGCBLin', 'Truth', 'Location', 'best');

% Decisions
% figure('Name', 'Decisions', 'NumberTitle', 'off');
subplot(2,2,2)
plot(1:t_final-1, bestSol, '-x');
hold on
plot(1:t_final-1, KGCBLinSol, '-o');
title('Decisions');
legend('MCTS', 'KGCBLin', 'Location', 'best');
xlabel('Iteration Counter');
ylabel('Decision');
axis([1 t_final-1 1 10])

% Accumulated Revenue
% figure('Name', 'Revenue', 'NumberTitle', 'off');
subplot(2,2,3)
plot(1:t_final-1, truth(bestSol), '-x');
hold on
plot(1:t_final-1, truth(KGCBLinSol), '-o');
title('Revenue');
legend('MCTS', 'KGCBLin', 'Location', 'best');
xlabel('Iteration Counter');
ylabel('Revenue');

% Cumulative Revenue
% figure('Name', 'Cumulative Revenue', 'NumberTitle', 'off');
subplot(2,2,4)
plot(1:t_final-1, cumsum(truth(bestSol)), '-x');
hold on
plot(1:t_final-1, cumsum(truth(KGCBLinSol)), '-o');
title('Cumulative Revenue');
legend('MCTS', 'KGCBLin', 'Location', 'best');
xlabel('Iteration Counter');
ylabel('Cumulative Revenue');
end