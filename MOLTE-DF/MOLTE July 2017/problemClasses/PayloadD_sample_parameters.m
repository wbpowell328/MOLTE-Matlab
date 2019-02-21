function theta = PayloadD_sample_parameters()
    min_E_ripe = 1.0;       % minimum value of E_ripe (energy barrier for diffusion/permeation)
    max_E_ripe = 1.18;       % maximum value of E_ripe
  

    K0_ripe = 10^11;
    K0_coal = 10^11;
    E_ripe = min_E_ripe + rand()*(max_E_ripe - min_E_ripe); 
    E_coal = E_ripe+ (1-2*rand())*0.1  ;

    theta(1) = K0_ripe;
    theta(2) = E_ripe;
    theta(3) = K0_coal;
    theta(4) = E_coal;
end


