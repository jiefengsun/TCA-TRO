function yarnBVP_ACT()
    % used for cable driven manipulator, do characterization. 
    %FE is the input force and young's modulus. 
    %Parameters
    % example. RodBVP_AB[2.4e3,0.02])
    E = 500e6; % 3.5e3; %Pa
    G = 400e6; % Pa 
    rho = 1.07e3; % kg/m^3 
    g = [0; 0; 9.81]; % 9.810 N/kg
    n_c = 10;
    r_c = 0.66e-3;
    L = n_c*2*pi*r_c; % m
    r = 0.35e-3;
    alpha = pi/8;
    rot_norm = 10; % rot/m
    A = pi*r^2;
    I = pi*r^4/4;
    J = 2*I;
    Kse = diag([G*A, G*A, E*A]);
    Kbt = diag([E*I, E*I, G*J]);
    N = 400; % Number of nodes
    v_0 = @(s) [0;0;1 ];

    %Boundary Conditions
    p0 = [0;0;0];
    R0 = eye(3);
    
    Me = [0;0;0]; % momentum at the tip, N*mm

    % Main Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    y = zeros(N,21);
    options = optimset('LargeScale','off','TolFun',1e-4,'MaxIter',25000,'MaxFunEvals',100000);

% %     init_guess =   y(1,13:18)';%% zeros(6,1); guess the initial m0 and n0, change the direction of the of the force, to get different solution. 
% 
%     fsolve(@RodShootingMethod, init_guess, options);%Use convex optimization to solve ICs 
%     end
    
   for ii =0:0.25:10
        rot_norm2 = rot_norm - ii;
        Fe =[0;0; 0]; %[0; 0; 0.01*i]; % force at the tip
        u_0 = @(s)  [0;0; 2*pi*rot_norm2/L;];
        
        init_guess = y(1,13:18)';%% zeros(6,1); guess the initial m0 and n0, change the direction of the of the force, to get different solution. 
                  
        fsolve(@RodShootingMethod, init_guess, options);%Use convex optimization to solve ICs
        
       % FM will return the initial values, which are the  correct guess for IC
       
        visual3D()
         Name = sprintf( 'fig%d.png', ii*4);
         Ndir = 'images\'; %'C:\MATLAB\RodDynamics_2020TRO\10_RodModel\06_TCA_yarn\images\'; 
         print(gcf, '-dpng','-r300',[Ndir, Name]) % this will create images in the current folder
         clf;
            
   end
    
 
    
    function visual3D()
        %Visualization
        fac = 1e-2;
        PG = zeros(N,3);
        PG2 = zeros(N,3);
        for i = 1:N
         PG(i,:) =  reshape(y(i,4:12),3,3)*[r*sin(pi/4);r*cos(pi/4);0]*0.5 + y(i, 1:3)';
         PG2(i,:) =  reshape(y(i,4:12),3,3)*[-r*sin(pi/4);-r*cos(pi/4);0]*0.5 + y(i, 1:3)';
        end
        y(:,1:3) = -y(:,1:3);
        PG = -PG;
        PG2 = -PG2; 
        plot3(y(:,1),y(:,2),y(:,3),'linewidth',5);  hold on
        plot3(PG(:,1),PG(:,2),PG(:,3),'r')
        plot3(PG2(:,1),PG2(:,2),PG2(:,3),'r')
        plot3(0,0,0,'b^', 'MarkerSize',6, 'LineWidth', 2)
        %plot3(y(:,19),y(:,20),y(:,21),'r');  hold on
        view(45, 30); 
        daspect([1 1 1]);

        %Farrow = [y(end, 1:3), y(end, 1:3) + fac*y(end,13:15)];
        %arrow3(Farrow(1:3), Farrow(4:6),'r')

        title('Untwisting of a Yarn');  axis([-r_c*5 r_c*5 -r_c*5 r_c*5 -0.6*L 0])
        grid off;    xlabel('x (m)');  ylabel('y (m)');  zlabel('z (m)');
        set(gca, 'color', 'none'); 
        % Enlarge figure to full screen.
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,0, .35, 0.8]);
    end

        
    %Subfunctions
    function residual = RodShootingMethod(guess) %Optimization objective function
        n0 = guess(1:3);                  %Update guessed initial conditions
        m0 = guess(4:6);
        
        y0 =  [p0; reshape(R0,9,1); n0; m0;[0;0;0]];

        
        [s,y] = ode45(@RodODE,linspace(0, L,N),y0);  %Numerically solve the resulting IVP
      
        nL_shot = y(end,13:15)';    % this is the position of the end for each tril.       
        mL_shot = y(end,16:18)';    % Calculate distal constraint violation
       
        Force_error = nL_shot - Fe;
        Moment_error = mL_shot - Me;
        residual = [Force_error ; Moment_error];
    end

    function ys = RodODE(s,y)
        R = reshape(y(4:12),3,3);
        n = y(13:15);
        m = y(16:18);

        v = Kse^-1*R.'*n + v_0(s);% [0;0;1]; Kirchoff rod
        u = Kbt^-1*R.'*m + u_0(s);

        ps = R*v;
        Rs = R*hat(u);
        ns = -rho*A*g;
        ms = -hat(ps)*n;
        PGs = Rs*[r;0;0] + ps;
        ys = [ps; reshape(Rs,9,1); ns; ms; PGs];
    end

    function skew_symmetric_matrix = hat(y)
        skew_symmetric_matrix = [  0   -y(3)  y(2) ;
                                  y(3)   0   -y(1) ;
                                 -y(2)  y(1)   0  ];
    end

end