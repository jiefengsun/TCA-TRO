% 12/10/2020 love full model with a hanging weight
% 12/15/2020 Modify to make the loading process to be the first. 
% 12/18/2020 further revise the constitutive law.
function  x_out = love_dynamics(volt, mw)
    cd C:\MATLAB\RodDynamics_2020TRO\05_Experiments_simulation\08_dynamics_with_load
    Name = sprintf( 'love_dyna_disp_%dV%dg.txt', volt, mw); 
    Name2 = sprintf( 'time_Pin_%dV%dg.txt', volt, mw); 
    dt = 0.01;
    [l_t, l_star, r_star, alpha_star, theta_bar_star, r_t, N_c, alpha_min] = TCA_geo(mw); 
    [t, T] = thermal_model(r_t, l_t, dt, Name2);
    [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(25, mw); 
    N = length(T);     
    g = 9.81; 
    alpha = alpha_star; 
    l = zeros(N, 1); 
    theta_bar = theta_bar_star;

    x_old = 0;
    x_cur = 0; 
    x_new = 0; 
    x_out = zeros(N, 1);
    F_old = mw/1000*g; 
      
    % solve initial shape 
    F_lv = @(alpha) F_old - GJ*cos(alpha_star)^2/r_star^2*(sin(alpha)- sin(alpha_star)) ...
          + EI*tan(alpha)*cos(alpha_star)^2/r_star^2*(cos(alpha)- cos(alpha_star));                
    alpha_l = fsolve(F_lv, alpha_star);       
    l_in = l_t * (sin(alpha_l)); 
    
    for i = 1:N
        fprintf('%d/%d \n', i, N); 
        F_old = fsolve(@func, F_old); %solve and update the F_old
        x_old = x_cur;
        x_cur = x_new; 
        x_out(i) = x_cur;
%         visualize(sol.y);
%         update the twist 
%         clf;
    end
    
    cd C:\MATLAB\RodDynamics_2020TRO\05_Experiments_simulation\08_dynamics_with_load
    writematrix([t, x_out], Name);
    plot(t, x_out)
    
    
    function res = func(F_in)
%         
        mw_var = F_in/g*1000;  
        [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(T(i), mw_var); 
        F_lv =@(alpha) F_in - GJ*cos(alpha_star)^2/r_star^2*(sin(alpha)- sin(alpha_star)) ...
            + EI*tan(alpha)*cos(alpha_star)^2/r_star^2*(cos(alpha)- cos(alpha_star));                
        alpha_l = fsolve(F_lv, alpha_l);       
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
        
        l  = l_l-D_l_h; 
        x_new = l_in - l ; 
        Fext = mw/1000*(g - (x_new- 2*x_cur+x_old)/dt^2);
        res = Fext - F_in; 
    end
   
end    