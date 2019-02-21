function [ mu_est, choices, ereward_sum, recommendedArm] = policyRun(mu_0,beta_W,covM,samples, policies, mu, tuned, alpha)

%INPUT;
%mu_0:= a column vector that contains the prior means for each decision
%covM:= M-by-M covariance matrix
%beta_W:= variance of the error from observations, for each decision
%samples: pre-generated sample realizations for each alternative N times(M-by-N matrix) 
%policies: all policies we want to compare (a cell with length K)
%mu:the expected value of each alternative
%alpha: tunable parameter
%tune: whether tune the parameter or not: if tune=0, use the default value

%OUTPUT:
%mu_est: estimated mean after N measurement for each policy
%choices: K-by-N matrix, each row contains alternatives picked at each iteration
%         for one policy
%reward_sum: a column vector of length K, accumulative reward collected by each policy in
%            N time steps
%recommendedArm: the best arm found by each policy after the measurement
%               budget exhausted.
[M,N]=size(samples);
mu_est=repmat(mu_0,1,length(policies));
ereward_sum=zeros(1,length(policies));
choices=zeros(length(policies), M);
recommendedArm=zeros(1,length(policies));
for i=1: length(policies)
    [mu_est(:,i), choices(i,:), recommendedArm(i)]= eval('caller', [char(policies(i)),'(mu_0,beta_W,covM,samples, alpha(i), tuned(i))']);
    %calculate the expected totol reward obtained by each policy
    ereward_sum(i)=sum(choices(i,:).*mu');
end

end

