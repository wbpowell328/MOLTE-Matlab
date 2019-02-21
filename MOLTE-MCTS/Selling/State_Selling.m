classdef State_Selling
    
    properties
        t;
        X_e;
        X_u;
        V_x = 0;    % Post decision state value
        action;     % Tracks the decision taken
        realization;% Tracks the realization
        Omega_e;
        Omega_u;
        children;
        index;
        N = 0;
        N_x = 0;
        actionChildren = containers.Map;
        pmf;
        inv;
        acc_profit = 0;
    end
    
    methods
        
        function this = State_Selling()
        end
        
        function obj = PreDecisionState(obj, X, index)
            obj.X_u = X;
            obj.index = index;
        end
        
        function obj = PostDecisionState(obj, Omega, index)
            obj.Omega_u = Omega;
            obj.index = index;
        end
    end
end