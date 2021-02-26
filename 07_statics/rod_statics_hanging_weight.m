function rod_statics_hanging_weight(mw, T_max)
    % input: mw the hanging weight, T_max is the max temperature; 
    % 
    % The final plot shows the displacement x with respect to the
    % temperature. 

    % Inputs
    T = (25:5:T_max)'; 
    mg = mw/1000*9.8; % mw is the weight in grams
    Fe_range = linspace(0, mg, floor(mg/0.05));
    
    % Geometry property
    [l_t, l_star, r_star, alpha_star, ~, N, alpha_min] = TCA_geo(mw);
    
    % Use 10 coils to accelerate the simulation.
    N_sim = 10; N_per_coil = 20; 
    Ns = N_per_coil*(N_sim); % Number of nodes   
    N_scale = N/N_sim; l_t = l_t/N_scale; 
    l_star = l_star/N_scale; 
    
    % Material property
    [EI, EA, GJ, GA, ~] = TCA_moduli_creeped(25, mw);% 
    K = [EI, EI, GJ, GA, GA, EA]'; 
    
    % Reference twists
    xi_0 = [0; cos(alpha_star)^2/r_star; sin(2*alpha_star)/(2*r_star);  0; 0; 1]; 
    %Boundary Conditions

    p0 = [r_star; 0; 0];
    R0 = [ -1,             0,            0;
            0,  sin(alpha_star), -cos(alpha_star);
            0,  -cos(alpha_star), -sin(alpha_star)];
    h0 = rotm2quat(R0)';
    Me = [0; 0; 0]; % momentum at the tip, N*mm
    Fe = [0; 0; 0]; % force at the tip 
   
    % Main Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialization
    solinit = bvpinit(linspace(0,l_t,Ns),@init);
    options = bvpset('Stats','off','RelTol',1e-5, 'NMax', 5000);
    N_step = length(T);
    l = ones(N_step, 1)*l_star;
    N_step1 = length(Fe_range);
    % iterate to the equilibrium with the after hanging the weigth. 
    for i  = 1:N_step1  % commont this part if we can load the data
       fprintf('preloading %d/%d \n', i, N_step1);
       Fe(3)=  - Fe_range(i);       
       sol = bvp5c(@static_ODE,@bc1,solinit,options);
       solinit = bvpinit(sol,[0 l_t]);
       visualize(sol.y);
       clf;
    end

    alpha = asin(norm(sol.y(1:3, end) - sol.y(1:3, 1))/l_t);
 
    for i = 1:N_step
        fprintf('%d/%d \n', i, N_step); 
         % update parameters 
        [EI, EA, GJ, GA, D_theta_bar_h] = TCA_moduli_creeped(T(i),  mw);% V_f = 0.4
        K =[EI, EI, GJ, GA, GA, EA]';
        % detect contact
       
        if(alpha < alpha_min)
             D_theta_bar_h = D_theta_bar_h +  20*exp(50*(alpha_min-alpha)); 
             fprintf('Reach the minimum length'); 
        end
    
        xi_0 =  [    0;
                          cos(alpha_star)^2/r_star;
                          sin(2*alpha_star)/(2*r_star)+ D_theta_bar_h ; 0;0;1];
        solinit = bvpinit(sol,[0 l_t]);
        sol = bvp5c(@static_ODE,@bc1,solinit,options); % used constrained BC
        
        % The solution at the mesh points
        ysol = sol.y; % we don't need arclength anymore. 
        l(i) = norm(ysol(1:3, end) - ysol(1:3, 1));
        alpha = asin(l(i)/l_t); %** this is IMPORTANT, update 
        visualize(ysol);
        % update the twist 
        clf;
    end
    x = -N_scale*(l-l(1));
    
   % writematrix([T, x], Name);
    plot(T,  x); 
    xlabel('Temperature $T$  ($^o$C)','interpreter','latex');
    ylabel('Displacement $x$  (mm)','interpreter','latex');
    
        
    function y = init(s)
        % initial shape
        y = [   r_star*cos((s*cos(alpha_star))/r_star)
                r_star*sin((s*cos(alpha_star))/r_star)
                s*sin(alpha_star)
                h0
                Me
                Fe];
    end
    
    function res = bc1(ya,yb)
         % boundary value
        res = [ ya(1:7) - [p0; h0]             
        Adg(yb(1:7))* yb(8:13) - [Me; Fe]];
    end
        

    function ys = static_ODE(s,y)
        % ODE     
        h = y(4:7);
        R = h2R(h);
        W = y(8:13); % W  = [m,n] First rotate 
        xi = K.^-1.*W + xi_0; 
        u = xi(1:3); 
        v = xi(4:6);
        
        ps = R*v;
        hs = h_diff(u, h);
        Ws = ad(xi)'*W; % ignore the distribued force
        ys = [ps; hs; Ws];        
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

  
    function visualize(y)
        % visualization
        plot3(y(1,:),y(2,:),y(3,:)); hold on
        plot3(y(1,end),y(2,end),y(3,end), 'ro', 'MarkerSize',10)
        title('TCA Dynamics');
        xlabel('x (m)');
        ylabel('y (m)');
        zlabel('z (m)');
        axis([-r_star*5 r_star*5 -r_star*5 r_star*5 -1.5*l_star 0]);
        grid on;
        % view(0,0)
        daspect([1 1 1]);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [.2,0.2, .25, .5]);
        drawnow;
    end
    
    function Adg = Adg(ph)
        % Adjoint representation of Lie group
        % used to transfer the 
        % ph is a vector ph = [p; h]; 
        p = ph(1:3); 
        h = ph(4:7);
        R = h2R(h);        
        Adg = [R, zeros(3,3); skew(p)* R, R];          
    end

end