function [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(T, mw)
% 10/05/2020 creat this script
% Cannot pass varialbe to dynamic work space (cannot load it using run), 
% Wrap it as a function and pass varialbes as arguments. - It seems
% clearer. 
% 10/28/2020 modified to in coporate the temperature change
% 12/19/2020 add the strees dependence.
% 1/2/2021 revided TCA_geo
 
    % import geometry
    [~, ~, ~, alpha_star, r_t_star, ~, ~] = TCA_geo(mw); 

    theta_bar_star = 2*pi*720/0.96; % rad/m
    
     % the polynomial for untwisting gamma     
    p_Gamma = [3.542e-06,-6.70546321037470e-05,1.00169272047785];
    Gamma = polyval(p_Gamma, T); % Gamma is related to the temperature. 
    r_t = r_t_star*Gamma;% consider the thermal expansion of the radius of the fiber. 
       
    % calculate the strain       
    G_ref = 0.22e9; 
    mg = mw*9.81/1000; 
    gamma_0 = 2*mg*cos(alpha_star)/(pi*r_t^3*G_ref); % strain    
    V_f = 0.35*0.887;
    p_gamma_0 = [3.24e-04,-0.027783,1.2]; %[134196665212.869,-487133894.722716,1104941.79219646];
    E_f = 3.9e9 - 1.1e6*(T-25)*polyval(p_gamma_0, gamma_0);    
    alpha_t = 0.8;
    E = 3*V_f*E_f/4*(1+cos(alpha_t))^2/(1+cos(alpha_t)+cos(alpha_t)^2);
    G = E_f*V_f/( pi*(1-cos(alpha_t))*sin(alpha_t)^3/(6*(alpha_t/2-1/4*sin(2*alpha_t))^2)  ...
             + (8*sin(alpha_t)^3)/(3*pi*(1-cos(alpha_t))*(1+cos(alpha_t))^2) ...
             + pi*(4-3*cos(alpha_t)-cos(alpha_t)^3)/(6*(alpha_t)/2-1/4*sin(2*alpha_t)*(1+cos(alpha_t)))); 
    A_t = pi*r_t^2; I = pi*r_t^4/4;  J = pi*r_t^4/2;
    EI = E*I;
    GJ = G*J;
    GA = G*A_t; % transverse shear
    EA = E*A_t;
        
   % untwisting parameter    
    D_theta_bar =  -theta_bar_star *(1-1./Gamma);. 
   
end