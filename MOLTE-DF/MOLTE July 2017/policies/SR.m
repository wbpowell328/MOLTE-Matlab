%{
notation for the following:
M is the number of alternatives.
N is the number of time-steps
M x N stands for a matrix with M rows and N columns

INPUT:
mu_0:   prior for the mean (M x 1), not used in this policy
beta_W: measurement precision (1/lambda(x)) (M x 1), not used in this policy
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

function [ mu_est, count,recommendedArm ] = SR( mu_0,beta_W,covM,samples,  alpha, tune)
%N>>M
[M,N]=size(samples);

count=zeros(M,1);  % count the times each alternative is measured

mu_est=zeros(M,1);

logBar=1/2+sum(1./(2:M));

surviving_arms=ones(M,1);

for m=1:M-1 % run for M-1 phases
  if (m==1)
      n_0=0;
  else
      n_0=n_1;
  end
  n_1=ceil(1/logBar*(N-M)/(M+1-m));
  count(surviving_arms==1)=count(surviving_arms==1)+n_1-n_0; %select each surviving arms for n_1-n_0 times
  observations=samples(surviving_arms==1, n_0+1:n_1);
  mu_est(surviving_arms==1)=(n_0*mu_est(surviving_arms==1)+sum(observations,2))/n_1;
  
  % remove one element from current surviving arms with the min mu_est, if there is a tie, select randomly 
  minindex=find(mu_est==min(mu_est(surviving_arms==1)) & surviving_arms==1);
  %randomly break the tie
  discardedArm=minindex(randi(length(minindex)));
  surviving_arms(discardedArm)=0;

end

recommendedArm=find(surviving_arms==1); %the recommended arm is the last surviving arm.

end
