function rod_dynamics_ode45(volt, mw)

    % This is the TCA static, 
    % the external force is prescribed as an boundary
    % condition. The force is increased step by step. 
    % using a bvp4c to solve the BVP problem
    % 10/05/2020 copied from function TCA_pass_stat_4c.m
    % borrow something from TCA_act_stat.m 
    % change all alpha_c to alpha
    % 10/6/2020 considering change of the coil radius
    % 10/7/2020 simulation works. Matched results. 
    % 10/28/2020 
    % 12/7 modified for hanging weight simulation
    % 12/8 include the fiber actuation model
    % 12/10 modify to make a general code for different weights. 
    % 12/14 roll back to the old constitutive law. 
    % 12/17/2020 modified for the rod varying load and incorporate the
    % thermal model. 
    % 12/20/2020 used new constitutive law. 
    % 1/5/2021 using ode45 to solve the dynamics
    Name = sprintf( 'rod_ode45_dyna_disp_%dV%dg.txt', volt, mw); 
    Name2 = sprintf( 'time_Pin_%dV%dg.txt', volt, mw); 
    
    mg = mw*9.81/1000; 
    F = [0:0.01:mg, mg]; % preloading 3g
    dt = 0.05;
    % independent of the fiber length 
    [l_t, l_star, r_star, alpha_star, theta_bar_star, r_t, N_c, alpha_min] = TCA_geo(mw); 
    [t, T] = thermal_model(r_t, l_t, 0.05, Name2);
  
    
    %   Name2 = sprintf( 'sol_%d.mat', mw); 
    N_sim = 10; % to scale the TCA down to accelerate the simulation. 
    N_per_coil = 20;
    Ns = N_per_coil*(N_sim); % Number of nodes   
    N_scale = N_c/N_sim; 
    l_t = l_t/N_scale; 
    l_star = l_star/N_scale; 
    % varying loading info
    [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(25, mw);% 
    alpha_star = - alpha_star; 
    %rho = 1.07e3; % kg/m^3 
    %g = 0*[0; 0; -9.81]; % 9.810 N/kg
    Kse = diag([GA, GA, EA]);
    Kbt = diag([EI, EI, GJ]);
    v_0 = [0; 0; 1];   
    u_0 = @(s)  [ (sin(theta_bar_star*s)*cos(alpha_star)^2)/r_star; 
                  (cos(theta_bar_star*s)*cos(alpha_star)^2)/r_star;
                  sin(2*alpha_star)/(2*r_star) ];
    %Boundary Conditions
    p0 = [r_star; 0; 0];
    R0 = [ -1,                0,               0;
            0, -sin(alpha_star), cos(alpha_star);
            0,  cos(alpha_star), sin(alpha_star)];
    h0 = rotm2quat(R0)';
    Me = [0; 0; 0]; % momentum at the tip, N*mm
    Fe = [0; 0; 0]; % force at the tip 
    theta_bar = theta_bar_star; 
    
    % Main Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    solinit = bvpinit(linspace(0,l_t,Ns),@init);
    options = bvpset('Stats','off','RelTol',1e-5, 'NMax', 5000);
    N_step = length(T);
    x_out = zeros(N_step, 1);
    alpha = -alpha_star; 
    N_it = length(F);

    for i  = 1:N_it % preloading
       Fe(3)=  -F(i);
       sol = bvp5c(@static_ODE,@bc1,solinit,options);
       solinit = bvpinit(sol,[0 l_t]);
%         visualize(sol.y);
%         clf;
    end
%     save(Name2,'-struct', 'sol')
%     sol = load(Name2); 



    x_old = 0;
    x_cur = 0; 
    x_new = 0; 
    F_old = -Fe(3);
    l_ref = norm(sol.y(1:3, end) - sol.y(1:3, 1));
    
    for i = 1:N_step
        fprintf('%d/%d \n', i, N_step); 
        F_old = fsolve(@func, F_old); %solve and update the F_old
        x_old = x_cur;
        x_cur = x_new; 
        x_out(i) = x_cur;
%         visualize(sol.y);
%        update the twist 
%         clf;
        
    end
    x_out = N_scale*x_out; 
    cd C:\MATLAB\RodDynamics_2020TRO\05_Experiments_simulation\08_dynamics_with_load
    writematrix([t, x_out], Name);
    plot(t, x_out)
    
    
    
    function res = func(F_in)
         % update parameters 
            Fe(3) = - F_in; 
            mw_var = F_in/9.8*1000;  
            [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(T(i), mw_var); 
            Kse = diag([GA, GA, EA]);
            Kbt = diag([EI, EI, GJ]);

            % detect contact
            if(alpha > alpha_min)
                 theta_bar = theta_bar_star +  D_theta_bar;
            else
                 theta_bar = theta_bar +  dt*10*exp(50*(alpha-alpha_min)); 
                 fprintf('Reach the minimum length'); 
            end

            u_0 = @(s)  [    (sin(theta_bar_star*s)*cos(alpha_star)^2)/r_star
                              (cos(theta_bar_star*s)*cos(alpha_star)^2)/r_star
                              sin(2*alpha_star)/(2*r_star)+ theta_bar ];
            solinit = bvpinit(sol,[0 l_t]);
            sol = bvp5c(@static_ODE,@bc1,solinit,options); % used constrained BC

            % The solution at the mesh points             
            x_new =  l_ref - norm(sol.y(1:3, end) - sol.y(1:3, 1));
            Fext = mw/1000*(9.81 - (x_new- 2*x_cur+x_old)/dt^2);
            res = Fext + Fe(3); 
    end
    
    function y = init(s)
        y = [   r_star*cos((s*cos(alpha_star))/r_star)
                r_star*sin((s*cos(alpha_star))/r_star)
                s*sin(alpha_star)
                h0
                Fe
                Me];
    end
    
    function res = bc1(ya,yb)
        res = [ ya(1:7) - [p0; h0] 
        yb(8:13) - [Fe; Me]];
    end
        

    function ys = static_ODE(s,y)
        % this is still required if the initial state is preloaded. 
        % This is another formulation using m and n as
        % varilables 
        h = y(4:7);
        R = h2R(h);
        n = y(8:10);
        m = y(11:13);

        v = Kse^-1*R.'*n + v_0;
        u = Kbt^-1*R.'*m + u_0(s);
        
        ps = R*v;
        hs = h_diff(u, h);
        ns = zeros(3,1); % ignore the body force
        ms = -cross(ps,n);

        ys = [ps; hs; ns; ms];
        
    end

    
    function R = h2R(h)
        %Quaternion to Rotation for integration
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
    end

    function hs = h_diff(u, h)
        % Calculate the derivative of the quaternion. 
        % if u to hs, if w we get ht
             hs = [ 0,    -u(1), -u(2), -u(3);
                    u(1),   0  ,  u(3), -u(2);
                    u(2), -u(3),   0  ,  u(1);
                    u(3),  u(2), -u(1),   0  ] * h/2;
    end

    %Function Definitions
    function visualize(y)
        plot3(y(1,:),y(2,:),y(3,:)); hold on
        plot3(y(1,end),y(2,end),y(3,end), 'ro', 'MarkerSize',10)
        title('TCA Dynamics');
        xlabel('x (m)');
        ylabel('y (m)');
        zlabel('z (m)');
        axis([-r_star*5 r_star*5 -r_star*5 r_star*5 -2.5*l_star 0]);
        grid on;
        % view(0,0)
        daspect([1 1 1]);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.45,0, .55, 1]);
        drawnow;
    end
  

end