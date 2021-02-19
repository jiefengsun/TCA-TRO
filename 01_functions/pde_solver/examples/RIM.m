function RIM
% These are Riemann problems for the Euler equations of gas dynamics 
% solved in H. Nessyahu and E. Tadmor, Non-oscillatory central differencing 
% for hyperbolic conservation laws, J. Comput. Phy. 87 (1990) 408-461. 

global gamma 
gamma = 1.4;

meth = menu('Specify the method','LxF','LxW','SLxW','NT');
switch meth
    case 1
        method = 'LxF';
    case 2
        method = 'LxW';
    case 3
        method = 'SLxW';
    case 4
        method = 'NT';
end

t = 0;
Npoints = 200;  
x = linspace(0,1,Npoints);

% Initial conditions
IC = menu('Specify initial data','Sod','Lax');
u = zeros(3,Npoints);
if IC == 1
    % Initial data of Sod
    problem = 'Sod';
    for m = 1:Npoints
        if x(m) < 0.5
            u(:,m) = [1,0,2.5]';
        else
            u(:,m) = [0.125,0,0.25]';
        end
    end
elseif IC == 2
    % Initial data of Lax
    problem = 'Lax';
    for m = 1:Npoints
        if x(m) < 0.5
            u(:,m) = [0.445,0.311,8.928]';
        else
            u(:,m) = [0.5,0,1.4275]';
        end
    end
end

NeumannL = 1:3;
NeumannR = 1:3;
sol = setup(3,@cl,t,x,u,method,false,[],{NeumannL,NeumannR});

% Plot density, velocity, pressure:
w = u;
w(2,:) = u(2,:)./u(1,:);                            % velocity 
w(3,:) = (gamma-1)*(u(3,:) - 0.5*u(2,:).*w(2,:));   % pressure
plot(x,w)
title(['Initial data for ',problem,' problem.'])
legend('density','velocity','pressure')
xlabel(['Solution at t = ',num2str(t),'.']);
if IC == 1
    axis([0 1 -0.5 2])
elseif IC == 2
    axis([0 1 -0.5 4])
end
pause

howfar = 0.04;
for m = 1:4
    sol = hpde(sol,howfar,@timestep);
    
    % Plot density, velocity, pressure:
    t = sol.t;
    u = sol.u;
    w = u;
    w(2,:) = u(2,:)./u(1,:);                            % velocity 
    w(3,:) = (gamma-1)*(u(3,:) - 0.5*u(2,:).*w(2,:));   % pressure
    plot(x,w)
    title([problem,' problem solved with ',method,'.'])
    legend('density','velocity','pressure')    
    xlabel(['Solution at t = ',num2str(t),'.']);
    if IC == 1
        axis([0 1 -0.5 2])
    elseif IC == 2
        axis([0 1 -0.5 4])
    end    
    pause

end

%=========================================================================
% Subfunctions
    
    % The variables are rho = u(1,:); m = u(2,:), E = u(3,:).  The variable
    % m = rho*v, where v is the velocity.  The NT paper writes the equations
    % as v_t + f(v)_x = 0, but hpde3 expects the form u_t = F(u)_x. 
    % Equation of state: p = (gamma-1)*(E - 0.5*rho*v^2)

    function F = cl(u)
        global gamma
        v = u(2,:)./u(1,:); 
        p = (gamma-1)*(u(3,:) - 0.5*u(2,:).*v);
        F = -[u(2,:); v.*u(2,:)+p; v.*(u(3,:)+p)];
    % end function cl
    
    function dt = timestep(dx,t,x,u)
        global gamma
        % sound speed is sqrt(gamma*p/rho).
        v = u(2,:)./u(1,:); 
        p = (gamma-1)*(u(3,:) - 0.5*u(2,:).*v);
        c = sqrt(gamma*p ./ u(1,:));
        dt = 0.9*dx/max( abs(v) + c );
    % end function timestep  

% end function RIM