function rod_conical_spiral_other_end()

    % This is the TCA static, 
    % the external force is prescribed as an boundary
    % condition. The force is increased step by step. 
    % using a bvp4c to solve the BVP problem
    % 10/05/2020 copied from function TCA_pass_stat_4c.m
    % borrow something from TCA_act_stat.m 
    % change all alpha_c to alpha
    % 10/6/2020 considering change of the coil radius
    % 10/7/2020 simulation works. Match results. 
    % 11/9/2020 copyied from rod_coil_conv.m
    % 1/26/2021 try to simulate conical spiral shape
    % 1/28/2021 simulate the case, the big coil end is fixed.  
    
    T = (25:5:100)'; 
    % add 06_TCA_Parameter in library in order to simulate. 
    [l_t, l_star, r_star, alpha_star, theta_bar_star, r_t, N_c, alpha_min] = TCA_geo(0); % input Number of coils and starting angle in Deg
    [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(25, 0); 
   
    
    a = 2*pi/0.005;
    b = tan(6/180*pi); % this is the b?
    rho = 1.07e3; % kg/m^3 
    g = 0*[0; 0; -9.81]; % 9.810 N/kg     
    Kse = diag([GA, GA, EA]);
    Kbt = diag([EI, EI, GJ]);
    Ns = 800; % Number of nodes    

    v_0 = [0; 0; 1]; 
    u_0 = @(z, s, theta_bar)  [0;
                               a*b*sqrt(4+a^2*z^2 + b^2*(2+a^2*z^2)^2)/(1+b^2*(1+a^2*z^2))^1.5;
                               a*(6+a^2*z^2)/(4+a^2*z^2+b^2*(2+a^2*z^2)^2) + theta_bar];
    %Boundary Conditions
    p0 = [0; 0; 0];
    R0 = [ -1, 0, 0;
            0, 0, 1;
            0, 1, 0];
    h0 = rotm2quat(R0)';
    Me = [0; 0; 0]; % momentum at the tip, N*mm
    Fe = [0; 0; 0]; % force at the tip 
    % Main Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname = [mfilename, '.avi']; aviObject = VideoWriter(fname);
    aviObject.FrameRate = 30; open(aviObject);
       
    solinit = bvpinit(linspace(0,l_t,Ns),@init);
    options = bvpset('Stats','off','RelTol',1e-5, 'NMax', 5000);
    N_step = length(T);
   
   for i = 1:N_step
        fprintf('%d/%d \n', i, N_step); 
        [EI, EA, GJ, GA, D_theta_bar] = TCA_moduli_creeped(T(i), 0); 
        Kse = diag([GA, GA, EA]);
        Kbt = diag([EI, EI, GJ]);
        theta_bar = theta_bar_star -  D_theta_bar; 
             
        sol = bvp5c(@static_ODE,@bc1,solinit,options);
        solinit = bvpinit(sol,[0 l_t]);
        % The solution at the mesh points
        ysol = sol.y; % we don't need arclength anymore. 
               
        visualize(ysol);
        % update the twist 
         
        for mm = 1:4; writeVideo(aviObject, getframe(gcf)); end; hold off;
        clf;
   end
  % writematrix([D_theta_range, l], 'rod_coil_conversion.txt');
    
   close(aviObject);    
   winopen( fname)
    
    function y = init(s)
        func =@(z) spiral_arclength(z) - s;
        z= fzero(func, 0);
        
        y = [   b*z*sin(a*z)
                b*z*cos(a*z)
                z
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
        
        func =@(z) spiral_arclength(z) - s;
        z = fzero(func, 0);
        v = Kse^-1*R.'*n + v_0;
        u = Kbt^-1*R.'*m + u_0(z,s,theta_bar); 
        
        ps = R*v;
        hs = h_diff(u, h);
        ns = zeros(3,1);
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
        axis([-0.01, 0.01, -0.05, 0.05, -0.01, 0.01]);
        grid on;
        % view(0,0)
        daspect([1 1 1]); % make axis equal
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.45,0, .55, 1]);
        drawnow;
    end



    function s = spiral_arclength(z) % arclength of conical spiral
       s = 1/2*z*sqrt(1+b^2*(1+a^2*z^2)) + (1+b^2)/(2*a*b)*asinh(a*b*z/sqrt(1+b^2));      
    end
  

end