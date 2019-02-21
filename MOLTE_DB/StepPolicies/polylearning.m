function [alpha] = polylearning(alphanought, n, beta)
% Polynomial learning rate rule 

% Input:
% alphanought - scaling estimate 
% n - iteration number 
% beta ? ? (1/2, 1] 

% Output 
% alpha - stepsize at iteration n (n != 0)

    if (n == 0) 
        error('N must be greater than zero');
    end 
    alpha = alphanought * 1 /(n.^beta);

end


