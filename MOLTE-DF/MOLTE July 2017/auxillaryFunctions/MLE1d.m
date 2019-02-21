%MLE for 1d case
function [alpha,beta,ave] = MLE1d( a , g )

    n=length(g);
    %define global variables for the use of patternsearch
    global X f;
    X=zeros(n);
    f=g;
    
    for i=1:n
        for j=1:n
            X(i,j)=(a(i)-a(j))^2;
        end
    end
    
    %find the parameter values that maxmize the likelihood
    a=patternsearch(@likelihood_1d, [10, 0.2],[],[],[],[],[1e-8,0]);

    alpha=a(1);
    beta=a(2);
   

    sigma=alpha*exp(-beta*X);
    ave=f'*(sigma\ones(n,1))/(ones(1,n)*(sigma\ones(n,1)));


    mean(f)
    
end

%calculate the negtive log-likelihood for one dimensional case
function w=likelihood_1d(parameter)
    global X f;
    alpha=parameter(1);
    beta=parameter(2);
        
    sigma=alpha*exp(-beta*X);
    n=length(f);

    
    ave=f'*(sigma\ones(n,1))/(ones(1,n)*(sigma\ones(n,1)));

    w=log(abs(det(sigma)))+(f-ave*ones(n,1))'*(sigma\(f-ave*ones(n,1)));
end



