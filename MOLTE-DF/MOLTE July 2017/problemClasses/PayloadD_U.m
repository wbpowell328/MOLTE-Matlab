function retVal = PayloadD_U(phi_w, phi_o, theta)
    % Constants
    T_f = 350;  % excited temperature
    T_0 = 290;  % normal tempreature
    tau_0 = 1;  % time scale under normal conditions
    tau_f = 1;  % time scale under excited conditions

    % evaluate utility function
    retVal = PayloadD_N(phi_w, phi_o, T_f, tau_f, theta)./PayloadD_N(phi_w, phi_o, T_0, tau_0, theta);
end



