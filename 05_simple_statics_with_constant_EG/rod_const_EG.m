function rod_const_EG()
    % This is the TCA static using constant E&G and ginores the conctat
    
    % The input
    mw = 30; % the weight hanged
    T = (25:5:160)'; % Temperature 
   
    % Geometry property     
    l_t =  175e-3; % Twistee fiber length
    r_start = 0.406e-3; % TCA diameter when it is made
    r_t = 0.22e-3; % twisted fiber radius
    alpha_star = 22.42/180*pi; % pitch angle
    n = l_t*cos(alpha_star)/(2*pi*r_start); % number of coils
    l_star = l_t*sin(alpha_star); %  TCA initial length
    r_star =  0.406e-3; % mm
  
    % Use 10 coils to accelerate the simulation.
    N_sim = 10; N_per_coil = 20; 
    Ns = N_per_coil*(N_sim); % Number of nodes   
    N_scale = n/N_sim; l_t = l_t/N_scale; 
    l_star = l_star/N_scale; 
    theta_bar_star = 4.71e+03; % rad/m
    
    % Material property
    E =  1.1980e+09; 
    G = 2.2177e+08;
    A_t = pi*r_t^2; I = pi*r_t^4/4;  J = pi*r_t^4/2;
    EI = E*I;   GJ = G*J;  GA = G*A_t;  EA = E*A_t;
    K = [EI, EI, GJ, GA, GA, EA]'; % using vector to speed up. 
    
    % Reference strains 
    xi_0 = [0; cos(alpha_star)^2/r_star; sin(2*alpha_star)/(2*r_star);  0; 0; 1]; 
    
    %Boundary Conditions
    p0 = [r_star; 0; 0];
    R0 = [ -1,             0,            0;
            0, sin(alpha_star), -cos(alpha_star);
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
    mg = mw/1000*9.8; % mw is the weight in grams
    Fe_range = linspace(0,mg, 30); 
    N_loading = length(Fe_range);
    % iterate to the equilibrium after hanging the weigth. 
    for i  = 1:N_loading   
       fprintf('Loading Step %d/%d \n', i, N_step); 
       Fe(3)=  -Fe_range(i);       
       sol = bvp5c(@static_ODE,@bc1,solinit,options);
       solinit = bvpinit(sol,[0 l_t]);
       visualize(sol.y);
       clf;
    end

    for i = 1:N_step
        fprintf('Actuation Step %d/%d \n', i, N_step); 
        % update parameters 
        p_Gamma = [3.54e-06,-6.7e-05,1]; %polynomial
        Gamma = polyval(p_Gamma, T(i)); % Gamma is related to the temperature. 
        % untwisting parameter    
        D_theta_bar_h = - theta_bar_star *(1-1./Gamma); 
        
        % update reference twist when heated
        xi_0 = [0;   cos(alpha_star)^2/r_star;
                sin(2*alpha_star)/(2*r_star)+ D_theta_bar_h;
                0;0;1];
        solinit = bvpinit(sol,[0 l_t]);
        sol = bvp5c(@static_ODE,@bc1,solinit,options); % used constrained BC
        
        % The solution at the mesh points
        ysol = sol.y; 
        l(i) = norm(ysol(1:3, end) - ysol(1:3, 1));
        visualize(ysol);
        clf;
    end
    x = N_scale*(l(1)-l);
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