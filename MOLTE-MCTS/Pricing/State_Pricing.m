classdef State_Pricing
    
    properties
        t;          % Time
        X_e;        % Explored decisions
        X_u;        % Unepxplored decisions
        V_x = 0;    % Post decision state value
        action;     % Tracks the decision taken
        realization;% Tracks the realization
        Omega_e;    % Explored observations
        Omega_u;    % Unexplored observations
        children;   % Children node
        index;      % Index in Tree
        N = 0;      % Number of times traversed (pre-decision state)
        N_x = 0;    % Number of times traversed (post-deicision state)
        actionChildren = containers.Map; % Links action to index of children in Tree
        theta;      % Coefficients
        B;      
        C;          % Covariance Matrix
        pmf;
        accR = 0;   % Accumulated revenue
    end
    
    methods
        
        function this = State_Pricing()
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