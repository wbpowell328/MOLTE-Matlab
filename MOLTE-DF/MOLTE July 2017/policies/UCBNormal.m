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


function [ mu_est, count, recommendedArm ] = UCBNormal( mu_0,beta_W,covM,samples,alpha, tune)
%N>>M
[M,N]=size(samples);

count=ones(M,1);  % count the times each alternative is measured
mu_est=samples(:,1);


V=zeros(M,1); % empirical variance
for i=M+1:N
    [minv, minindex]=min(count);
    if (minv < ceil(8*log(i)))  % if an arm is played less than ceil(8log(i)),
         x=minindex;            % choose to play this arm
    else %else choose the arm with the largest UCB index        
      UCBindex=mu_est+sqrt(16*V).*sqrt(log(i-1)./count);
      [max_est, x]=max( UCBindex);
    end
    count(x)=count(x)+1;
    W=samples(x, count(x));
    temp=mu_est(x);
    mu_est(x)=((count(x)-1)*mu_est(x)+W)/count(x);
    V(x)=(1-1/(count(x)-1))*V(x)+count(x)*(mu_est(x)-temp)^2;

end
[value, recommendedArm] = max(mu_est);
end

