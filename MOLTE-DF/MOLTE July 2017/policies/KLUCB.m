%{
notation for the following:
M is the number of alternatives.
N is the number of time-steps
M x N stands for a matrix with M rows and N columns

INPUT:
mu_0:   prior for the mean (M x 1), not used in this policy
beta_W: measurement precision (1/lambda(x)) (M x 1)
covM:   initial covariance matrix (M,M), not used in this policy
samples: pre-generated observation realizations (M x N)
alpha: tunable parameter
tune: whether tune the parameter or not: if tune=0, use the default value

OUTPUT:
mu_est:     Final estimates for the means (M x 1)
count:      The number of each alternative being measured (M x 1)
recommendedArm:  The index of the arm recommended by this policy after the
measurement budget exhausted
%}


function [ mu_est, count, recommendedArm] = KLUCB( mu_0,beta_W,covM,samples,alpha, tune)
%N>>M
[M,N]=size(samples);

count=ones(M,1);  % count the times each alternative is measured
mu_est=samples(:,1); % start by measuring each alternative once
rewards=mu_est.*ones(M,1); % initialize the reward for each arm

if tune==0 %not tuned
    alpha=0;
end


for i=M+1:N
    % calculate the upper confidence bound for each arm
    b=log(i)+alpha*log(log(i));
    q=sqrt(1./beta_W).*sqrt(2*b./count)+rewards./count;
    [max_est, x]=max(q);
    count(x)=count(x)+1;
    W=samples(x, count(x));
    rewards(x)=rewards(x)+W;
    mu_est(x)=((count(x)-1)*mu_est(x)+W)/count(x);
end

[value, recommendedArm] = max(mu_est);
end