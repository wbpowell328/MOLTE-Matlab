%% Pricing Problem
clear;

hold off
%
X_i1=[1;2;3;4;5;6;7;8;9;10];
X_i2=[1;4;9;16;25;36;49;64;81;100];

X=[ones(10,1),X_i1,X_i2];

%Prior Estimates
%theta_0=[64;2;-.2];
theta_0=[64;4;-.4];

%Standard deviations for Theta
%theta_Sigma=[20;5;.5];
theta_Sigma=[10;1;.1];

%Pick a true value for theta
thetaT=theta_0+randn(3,1).*theta_Sigma

%Set the true values as X*theta
truth=X*thetaT;

%Set the measurement variance
MVar=4^2;

%Run the KGCBLinR L many times
L=5;

%Set the prior on covariance matrix, assume independence
C=diag(theta_Sigma).^2;

% Parameters
t_final = 7;
d_thr = 10;
e_thr = 5;
alpha = 1;

stepsize = 0.01;

% weight between 0 - 1, higher weight gives bias on info_gain

plot(X_i1, truth);

% Initialize the B, C, and theta
num_init = 2;
Xs = X(randperm(length(X_i1),3),:);
Y = Xs*thetaT + normrnd(0, sqrt(MVar), 3, 1);

B = inv(Xs'*Xs);
theta = B*Xs'*Y;
solution = zeros(t_final-1,1);
for t = 1:t_final-1
    fprintf('t = %d: \n', t);
    clear Tree;
    clear childToParent;
    
    Tree(1)         = State_Pricing();
    Tree(1)         = Tree(1).PreDecisionState(1:10, 1);
    Tree(1).B       = B;
    Tree(1).theta   = theta;
    Tree(1).C       = C;
    Tree(1).t       = t;
    index = 1;
    childToParent(1) = 0;
    
    
    for i = 1:200
        [Tree, childToParent, index] = TreePolicy_Pricing(Tree, t_final, d_thr, e_thr, alpha, childToParent, thetaT);
        Tree(index) = SimPolicy_Pricing('abc', Tree(index), Tree(childToParent(index)).theta,  X, MVar, truth, L, stepsize, .6);
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
%         fprintf('Action: %d, Value: %d, Count: %d\n', currNode.action, currNode.V_x, currNode.N_x);
        
        %     hold on
        %     plot(X_i1, X*Tree(currNode.children(1)).theta);
        %     pause
    end
    fprintf('Action: %d, Value: %d, Count: %d\n', Tree(maxNode).action, Tree(maxNode).V_x, Tree(maxNode).N_x);
    theta = Tree(maxNode).theta;
    B = Tree(maxNode).B;
    C = Tree(maxNode).C;
    
    
end


%% Comparing to KGCBLin
[thetaLast,OC,choices,thetaEst]=KGCBLinR(X,theta,MVar,truth,t_final,C,0,1);