function[alpha, beta, v, lambda] = BAKF(iter, x, xprev, nu, beta, v, lambda, a)
    
    [beta,v,var,lambda] = updateStats(x, xprev, nu, beta, v, lambda, a);

    % update alpha for next time step
    if iter == 1 || iter == 2 
        a = (1) / (iter+1);
        alpha = a;
    end 
    
%     a = 1 - ((var) / ((1+lambda)*var + beta*beta));
    a = 1 - (var)/(v);
    alpha = a;

end 

function [beta1,v1,var1,lambda1] = updateStats(theta, prevtheta, nu, beta, v, lambda, a)
   
    epsilon = prevtheta - theta;
    beta1 = (1-nu)*beta + nu*epsilon;
    v1 = (1-nu)*v + nu*epsilon*epsilon;
    var1 = (v-beta*beta) / (1+lambda);
    lambda1 = (1 - a)*lambda + (a)^2;
    
end 