% 12/10/2020 cst full model with a hanging weight
% update the temperature and stress dependent modulus. 
% 12/15/2020 Modify the
% 2/12/2021 use this version and modified for uploading
function x = cst_hang_weight(mw)

   % Fe here is positive 
    Fe = mw/1000*9.8; 
    Name = sprintf( 'cst_free_stroke_%dg.txt', mw); 
    T = (25:5:160)'; 
    N_step = length(T); 
    [l_t, l_star, r_star, alpha_star,  r_t, N_c, alpha_min] = TCA_geo(mw); % input Number of coils and starting angle in Deg
    alpha = alpha_star; 
    alpha_h = zeros(N_step, 1); 
    D_l = zeros(N_step, 1);
    for i = 1:N_step
        [EI, EA, GJ, GA, D_theta_bar_h] = TCA_moduli_creeped(T(i), mw);
        
        if(alpha < alpha_min)
             D_theta_bar_h =  D_theta_bar_h +  10*exp(50*(alpha_min-alpha)); 
             fprintf('Reach the minimum length'); 
        end
       % This is to solve the 
        func1 = @(alpha)   l_t*(sin(alpha_star)- sin(alpha)) + l_t*r_star*cos(alpha)^2/cos(alpha_star)* D_theta_bar_h; 
        alpha_h(i) = fsolve(func1, alpha_star);
        
        
        func2 =@(alpha)  l_t*(sin(alpha)- sin(alpha_h(i))) - K_cst_inv(alpha)*Fe;
        alpha = fzero(func2,alpha_star );
        
        D_l(i) =  l_t*(sin(alpha_star)-sin(alpha)); 
 
    end
    x = D_l-D_l(1);
    writematrix([T, x], Name);
     
    
    
    function   K_cst_inv  = K_cst_inv(alpha)
       %
       r = r_star*cos(alpha)/cos(alpha_star);
       f_11 = l_t*(cos(alpha)^2/GA + sin(alpha)^2/EA + (r^2*sin(alpha)^2)/EI + (r^2*cos(alpha)^2)/GJ);
       K_cst_inv = f_11;       
    end    
    
    
end
    