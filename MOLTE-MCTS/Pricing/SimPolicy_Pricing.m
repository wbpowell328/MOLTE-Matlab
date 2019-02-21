function Node = SimPolicy_Pricing(rolloutPolicy, Node, prevTheta, X, MVar, truth, L, stepsize, weight, param)
% Runs the rollout policy and returns a value. In this problem set the
% value of a node as the difference between the two curves
% Inputs
%   rolloutPolicy   = function handle of rollout policy
%   Node            = the node which we are returning the value for
%   X               = the design matrix (the experiments)
%   MVar            = variance of measurement noise
%   truth           = truth
%   L               = number of measurements
% Output
%   Node            = node with an update value

if Node.V_x ~= 0
    return;
end

if L == 0
    disp(L)
    Node.V_x = (1-weight)*Node.accR;
    return;
end

[theta,choices] = rolloutPolicy(X,Node.theta,MVar,truth,L,Node.C,0,param);
choices = choices';

[m, n] = size(X);   %m: time periods, n: parameters
max_x = max(X(:,2)); % the x's, first column are 1s, 2..i..n columns are x^(i-1)
x = (0:stepsize:max_x)';
X_fine = ones(length(x), n);
X_choice = ones(length(choices), n);

for i = 2:n
    X_fine(:,i) = x.^(i-1);
    X_choice(:,i) = choices.^(i-1);
end

info_gain = norm(X_fine*abs(theta(:,end) - prevTheta))/length(X_fine); % average difference between curves (2-norm)

revenue = sum(X_choice*theta(:,end))/Node.t; % Approximation of average revenue 

Node.V_x = weight*info_gain + (1-weight)*(revenue + Node.accR);

end