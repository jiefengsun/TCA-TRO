function frame_demo_fd()
   
    % To demonstrate the frame transformation, using a finite differce solver. 
    

    % Geometry property     
    l_t =  175e-3; % Twistee fiber length
    r_start = 0.406e-3; % TCA diameter when it is made
    r_t = 0.22e-3; % twisted fiber radius
    alpha_star = 18/180*pi; % pitch angle
    n = l_t*cos(alpha_star)/(2*pi*r_start); % number of coils
    r_star =  0.406e-3; % mm
    h = 0; 
    % Use 10 coils to accelerate the simulation.
    N_sim = 2; N_per_coil = 50; 
    Ns = N_per_coil*(N_sim); % Number of nodes   
    N_scale = n/N_sim; l_t = l_t/N_scale; 
   
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
    arrow_color = {'b' 'g' 'r' 'y' 'c'};    %color of the basis vectors
    N_f = 1; 
    currentFolder = pwd;
    Ndir = [currentFolder,'\images'];
    if ~exist(Ndir, 'dir')
       mkdir(Ndir)
    end
    
    % Main Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialization
    solinit = bvpinit(linspace(0,l_t,Ns),@init);
    options = bvpset('Stats','off','RelTol',1e-5, 'NMax', 5000);

 
    % iterate to the equilibrium after hanging the weigth. 
    sol = bvp5c(@static_ODE,@bc1,solinit,options);
   
    plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:))
    axis([-1.5*r_star 1.5*r_star -1.5*r_star 1.5*r_star -5*r_star r_star]);
    for j= 1:4:Ns

       hold on
      frame =  plot_bf(h2R(sol.y(4:7, j)), sol.y(1:3,j));
      delete(frame);

    end
    
    
    function h =  plot_bf(R_b, p) % plot body frame
    
    %
    hold on
   % axis([-r_star r_star -r_star r_star -5*r_star r_star]);
    for i = 1:3     
      h(i)  = quiver3( p(1), p(2), p(3), R_b(1,i),R_b(2,i),R_b(3,i), 0.0003,'linewidth',2, 'color', arrow_color{i}); 
      h(i).MaxHeadSize = 30; 
      h(i).MarkerSize =  25;       
    end
    grid on
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1,0.1, .3, 0.5]);
    set(gca, 'color', 'none')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    daspect([1 1 1]);
    drawnow 
    Name = sprintf( '/fig%d.png', N_f);  print(gcf, '-dpng','-r200',[Ndir, Name])  
    N_f = N_f+1; 
   
    end

    
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
      Adg(yb(1:7))* yb(8:13) -  [Me; Fe]];
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
        Ws = ad(xi)'*W; % ignore the distribued wrench
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