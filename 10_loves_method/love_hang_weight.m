% 12/10/2020 love full model with a hanging weight
% 12/15/2020 Modify to make the loading process to be the first.  
function x = love_hang_weight(mw)
    
    Fe = mw/1000*9.8; 
    Name = sprintf( 'love_free_stroke_%dg.txt', mw); 
     
    T = (25:5:160)'; 
    N = length(T); 
    [l_t, l_star, r_star, alpha_star, theta_bar_star, r_t, N_c, alpha_min] = TCA_geo( mw); % input Number of coils and starting angle in Deg
    
    alpha = alpha_star; 
    l = zeros(N, 1); 
    alpha_l =  alpha_star; 
    for i = 1:N
        [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(T(i), mw); 
        F_lv =@(alpha) Fe - GJ*cos(alpha_star)^2/r_star^2*(sin(alpha)- sin(alpha_star)) ...
           + EI*tan(alpha)*cos(alpha_star)^2/r_star^2*(cos(alpha)- cos(alpha_star)); 
            
%         r_ratio =@(alpha) (sin(alpha_star)*cos(alpha_star)*tan(alpha) + EI/GJ*cos(alpha_star)^2)/(cos(alpha)^2*(EI/GJ + tan(alpha)^2)); 
%         F_lv = @(alpha) Fe - EI/r_star^2* r_ratio(alpha)*(cos(alpha_star)^2-r_ratio(alpha)*cos(alpha)^2)/sin(alpha); 
       alpha_l = fsolve(F_lv,  alpha_l );   
                
        l_l = l_t * (sin(alpha_l));
        
        
        if(alpha > alpha_min)
             theta_bar = theta_bar_star +  D_theta_bar ;            
        else
             theta_bar = theta_bar +  10*exp(50*(alpha-alpha_min)); 
             fprintf('Reach the minimum length'); 
        end
        D_theta = theta_bar*l_t; 
        A_lv = r_star/cos(alpha_star); 
        D_l_h = A_lv*D_theta; 
        
        l(i) = l_l-D_l_h; 
        alpha = asin(l(i)/l_t); % update alpha for the contact dectction
     
    end
    x = -(l-l(1)); 
    cd C:\MATLAB\RodDynamics_2020TRO\05_Experiments_simulation\02_Active_with_load\love_hang_sim
    writematrix([T, x ], Name);
    
        
end
    