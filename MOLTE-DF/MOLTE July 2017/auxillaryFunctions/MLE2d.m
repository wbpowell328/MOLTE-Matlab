%MLE for 2d case
function [alpha,beta1,beta2,ave] = MLE2d( a,b,g )

    n=length(g);
    %define global variables for the use of patternsearch
    global X Y f;
    X=zeros(n);
    Y=zeros(n);
    f=g;
    
    for i=1:n
        for j=1:n
            X(i,j)=(a(i)-a(j))^2;
            Y(i,j)=(b(i)-b(j))^2;
        end
    end
    
    %find the parameter values that maxmize the likelihood
    aa=patternsearch(@likelihood_2d,[10,0.2,0.2],[],[],[],[],[0.001, 0.0001,0.0001]);

    alpha=aa(1);
    beta1=aa(2);
    beta2=aa(3);
     
    sigma=alpha*exp(-beta1*X-beta2*Y);
    ave=f'*(sigma\ones(n,1))/(ones(1,n)*(sigma\ones(n,1)));


    mean(f);
end

%calculate the negtive log-likelihood for one dimensional case
function w=likelihood_2d(parameter)
    global X Y f;
    alpha=parameter(1);
    beta1=parameter(2);
    beta2=parameter(3);   
    
    sigma=alpha*exp(-beta1*X-beta2*Y);
    n=length(f);
    ave=f'*(sigma\ones(n,1))/(ones(1,n)*(sigma\ones(n,1)));

    w=log(abs(det(sigma)))+(f-ave*ones(n,1))'*(sigma\(f-ave*ones(n,1)));

end



