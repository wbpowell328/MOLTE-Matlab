%{
notation for the following:
M is the number of alternatives.
N is the number of time-steps
M x N stands for a matrix with M rows and N columns

INPUT:
mu_0:   prior for the mean (M x 1)
beta_W: measurement precision (1/lambda(x)) (M x 1)
covM:   initial covariance matrix (M,M)
samples: pre-generated observation realizations (M x N)
alpha: tunable parameter
tune: whether tune the parameter or not: if tune=0, use the default value

OUTPUT:
mu_est:     Final estimates for the means (M x 1)
count:      The number of each alternative being measured (M x 1)
recommendedArm:  The index of the arm recommended by this policy after the
measurement budget exhausted
%}


function [ mu_est, count, recommendedArm ] = ThompsonSampling( mu_0,beta_W,covM,samples, alpha, tune)

[M,N]=size(samples);

count=zeros(M,1);  % count the times each alternative is measured
mu_est=mu_0;

for i=1:N
    posteriorSample=mvnrnd(mu_est,covM);  % sample from the posterior distribution
    [max_est, x]=max(posteriorSample);
    count(x)=count(x)+1;
    W=samples(x, count(x));
    
    e_x=zeros(M,1);
    e_x(x)=1;
    
    %updating equations for Normal-Normal model with covariance
    addscalar = (W - mu_est(x))/(1/beta_W(x) + covM(x,x));
    mu_est=mu_est + addscalar*covM*e_x;
    covM = covM - (covM*(e_x*e_x')*covM)/((1/beta_W(x)) + covM(x,x));

end
[value, recommendedArm] = max(mu_est);
end
