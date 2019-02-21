% This function returns the amount payload delivered as a function of initial internal water
% concentration, oil droplet diameter, temperature and time (tau). The kinetic coefficients
% are held in theta.

function retVal = PayloadD_N(phi0_w, phi_o, T, tau, theta)

% unknown variables
K0_ripe = theta(1);
E_ripe = theta(2);
K0_coal = theta(3);
E_coal = theta(4);

% fixed constants
d_o = 9;
DT = 0.0001;        %time step to integrate over
N0 = 100;
V_o = pi/6.*d_o.^3;
V_e = V_o./phi_o;
k = 8.6173324e-5;


% initial conditions
V_i = phi0_w.*V_o;
N_i = N0;
phi_w = phi0_w;

% integrate
for t = 1:DT:tau
   K_ripe = pi.*K0_ripe.*d_o.^3.*phi_w.*(1-phi_o)./((1-phi_o) + phi_w.*phi_o).*exp(-E_ripe./(k.*T));
   K_coal = K0_coal.*phi_w.*phi_o.*exp(-E_coal./(k.*T));

   % Update
   N_i_new = N_i - K_coal.*N_i.*DT - K_ripe.*(N_i./V_i - (N0-N_i)./V_e)*DT;
   V_i_new = V_i - K_coal.*V_i.*DT;
   phi_w_new = V_i ./V_o;

   % Step
   N_i = N_i_new;
   V_i = V_i_new;
   phi_w = phi_w_new;

   % In case we overstepped 0 because of large DT:
   N_i = N_i.*(N_i > 0);
   V_i = V_i.*(V_i > 0);
   phi_w = phi_w.*(phi_w > 0);

end

% return amt payload in external phase
retVal = N0 - N_i;

end

