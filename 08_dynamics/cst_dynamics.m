% 12/10/2020 cst full model with a hanging weight
% update the temperature and stress dependent modulus. 
% 12/15/2020 Revise the code to first consider the loading process
function  x_out = cst_dynamics(volt, mw)

   % Fe here is positive      
    Name = sprintf( 'cst_dyna_disp_%dV%dg.txt', volt, mw); 
    Name2 = sprintf( 'time_Pin_%dV%dg.txt', volt, mw); 
    
    dt = 0.02;
    [l_t, l_star, r_star, alpha_star, theta_bar_star, r_t, N_c, alpha_min] = TCA_geo(mw); 
    [t, T] = thermal_model(r_t, l_t, dt,Name2);
    [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(25, mw); 
    N = length(T);
    g = 9.81; 
    alpha= alpha_star; 
    l = zeros(N, 1);
    theta_bar = theta_bar_star;
    x_old = 0;
    x_cur = 0; 
    x_new = 0; 
    x_out = zeros(N, 1);
    F_old = mw/1000*9.8; 
    
    func_l =@(alpha)  l_t*(sin(alpha)- sin(alpha_star)) - K_cst_inv(alpha)*(F_old);
    alpha_l = fzero(func_l,alpha_star );
    l_in = l_t * (sin(alpha_l)); 
    
    
    for i = 1:N
        fprintf('%d/%d \n', i, N); 
        F_old = fsolve(@func, F_old); %solve and update the F_old
        x_old = x_cur;
        x_cur = x_new; 
        x_out(i) = x_cur;
%       visualize(sol.y);
%       update the twist 
%       clf;
    end
    
    cd C:\MATLAB\RodDynamics_2020TRO\05_Experiments_simulation\08_dynamics_with_load
    writematrix([t, x_out], Name);
    plot(t, x_out)
    
    function res = func(F_in)
         mw_var = F_in/g*1000;  
        [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(T(i), mw_var); 
        
        func_l =@(alpha)  l_t*(sin(alpha)- sin(alpha_star)) - K_cst_inv(alpha)*(F_in);
        alpha_l = fzero(func_l,alpha_star );
        if(alpha > alpha_min)
             theta_bar = theta_bar_star +  D_theta_bar ;
        else
             theta_bar = theta_bar +  5*dt*exp(50*(alpha-alpha_min)); 
             fprintf('Reach the minimum length'); 
        end
        
        D_theta = theta_bar*l_t; 
        func_h = @(alpha)   l_t*(sin(alpha_l)- sin(alpha)) - r_star*cos(alpha)^2/cos(alpha_star)* D_theta; % This is only consider the f_11, case
        alpha = fsolve(func_h, alpha);
        
        l =  l_t* sin(alpha); 
        x_new = l_in - l ; 
        Fext = mw/1000*(g - (x_new- 2*x_cur + x_old)/dt^2);
        res = Fext - F_in;  
     end
    
    function   K_cst_inv  = K_cst_inv(alpha)
       % calculate the conpliant coefficient, which is the inverse of the
       % stiffness. 
       r = r_star*cos(alpha)/cos(alpha_star);
       f_11 = l_t*(cos(alpha)^2/GA + sin(alpha)^2/EA + (r^2*sin(alpha)^2)/EI + (r^2*cos(alpha)^2)/GJ);
       K_cst_inv = f_11;       
    end

    
end
    