function Node = SimPolicy_Selling(Node, t_final, forecast, var_ub, rolloutPolicy)

if Node.V_x ~= 0
    return;
end

D = forecast.demand;
D = round(D(Node.t:t_final));
P = forecast.price(Node.t:t_final);

if Node.N == 0 && Node.t <= t_final
    
    fval = rolloutPolicy(Node.inv, D, P, var_ub);

    % Assign Value to Node
    Node.V_x = Node.acc_profit + fval;
    
end
end