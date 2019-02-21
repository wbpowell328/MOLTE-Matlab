%% Rollout policy
% r       = starting inventory
% D       = demand
% P       = cost
% var_ub  = max production

    function [fval, decisions] = DeterministicIPSelling(r, D, P, var_ub)
        
        n_x     = length(D);    % Time periods forward
        D       = min(D, r);
        % Objective function
        % max Pt*(Dt - St) -> max Pt*Dt - Pt*St -> max -Pt*St (+ Constant)
        f       = zeros(1, 2*n_x);
        f(1:n_x) = 0;
        f(n_x+1:2*n_x) = -P;
        
        % RHS
        beq     = min(D, r);
        beq(1)  = beq(1) - r;
        
        % Constraint matrix
        e       = ones(n_x, 1);
        I       = spdiags([e -1*e], -1:0, n_x, n_x);
        uS      = spdiags(e, 0, n_x, n_x);
        Aeq     = [I uS];

        % Bounds
        vars    = ones(2*n_x,1);
        lb      = 0*vars;
        ub      = [var_ub*e; D];
        
        %% INTLINPROG IMPLEMENTATION
        %         options = optimoptions('intlinprog');
        %         options.Display = 'off';
        %         [x, fval, exitflag] = intlinprog(f, 1:n_x, [], [], Aeq, beq, lb, ub, options);
        
        %% GUROBI IMPLEMENTATION - works as same as intlinprog in MATLAB, aside from options
        [x, fval, exitflag] = intlinprog_GUROBI(-f, 1:n_x, [], [], Aeq, beq, lb, ub);
        
        decisions = D - x(n_x+1:end);
        
        fval = P*decisions;
        
    end