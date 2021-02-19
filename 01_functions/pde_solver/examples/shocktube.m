function shocktube
% HW 0.0.3 (Shock Tube Problem) of J.W. Thomas, Numerical Partial
% Differential Equations Finite Difference Methods, Springer, New York, 
% 1995. This is a system with a shock, contact discontinuity, and expansion 
% fan. Thomas discusses boundary conditions at length and then solves the
% problem with homogeneous Neumann boundary conditions.  

global gamma 

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
if meth == 4
    form = 3;
else
    form = menu('Specify the form of the PDEs','1','2','3');
end
switch form
    case 1
        pdefun = @pdes;
    case 2
        source = [0; 0; 0];
        pdefun = {@flux,source};
    case 3
        pdefun = @cl;
end

t = 0;
% Mesh as in plots of Thomas
Npoints = 200;    
x = linspace(-2,2,Npoints);

% Initial conditions
v0 = zeros(size(x));
rho0 = ones(size(x));
rho0(find(x < 0)) = 2;
% Thomas specifies the initial pressure and computes the initial 
% energy from the equation of state.
gamma = 1.4;
p0 = rho0;
E0 = p0/(gamma-1) + 0.5*rho0.*v0.^2;   
u = [rho0; rho0.*v0; E0];

NeumannL = 1:3;
NeumannR = 1:3;
sol = setup(form,pdefun,t,x,u,method,[],[],{NeumannL,NeumannR});

% A fixed time step as in plots of Thomas
timestep = 0.0025;

w = u;
w(2,:) = u(2,:)./u(1,:);                            % velocity 
w(3,:) = (gamma-1)*(u(3,:) - 0.5*u(2,:).*w(2,:));   % pressure
close all
plot(x,w)
title('Initial data')
legend('density','velocity','pressure')
xlabel(['Solution at t = ',num2str(t),'.']);
axis([-2 2 -0.2 2.2])
legend('density','velocity','pressure')
pause


howfar = 0.2;
for m = 1:5
    sol = hpde(sol,howfar,timestep);    
    % Plot density, velocity, pressure:
    t = sol.t;
    u = sol.u;
    w = u;
    w(2,:) = u(2,:)./u(1,:);                            % velocity 
    w(3,:) = (gamma-1)*(u(3,:) - 0.5*u(2,:).*w(2,:));   % pressure
    plot(x,w)
    title(['Computed with ',method,' and form ',num2str(form),'.'])
    legend('density','velocity','pressure')
    xlabel(['Solution at t = ',num2str(t),'.']);
    axis([-2 2 -0.2 2.2])
    pause
end
 
%=========================================================================
% Subfunctions

% Use formulas in Thomas, pp. 348-349.   In his notation v = [rho; m; E]. 
% Here rho = u(1,:); m = u(2,:), E = u(3,:).  The signs are changed in
% these functions because his formulas are written for v_t + f(v)_x = 0 
% and v_t + A*v_x = 0. Equation of state: p = (gamma-1)*(E - 0.5*rho*v^2)

    function F = cl(u)
        global gamma    
        v = u(2,:)./u(1,:); 
        p = (gamma-1)*(u(3,:) - 0.5*u(2,:).*v);
        F = -[u(2,:); v.*u(2,:)+p; v.*(u(3,:)+p)];
    % end function cl

    function F = flux(t,x,u)
        F = feval(@cl,u);
    % end function flux  

    function F = pdes(t,x,u,u_x)
        global gamma
        F = zeros(size(u));
        F(1,:) = u_x(2,:);
        F(2,:) = 0.5*(gamma-3)*(u(2,:).^2 ./ u(1,:).^2) .* u_x(1,:) ...
                -(gamma-3)*(u(2,:)./u(1,:)).*u_x(2,:) + (gamma-1)*u_x(3,:);
        F(3,:) = -(u(2,:)./u(1,:).^2).*(gamma*u(3,:) - (gamma-1)*...
                 (u(2,:).^2./u(1,:))).*u_x(1,:) + ((gamma*u(3,:) - ...
                 1.5*(gamma-1)*(u(2,:).^2./u(1,:)))./u(1,:)).*u_x(2,:)...
                 + gamma*(u(2,:)./u(1,:)).*u_x(3,:);
        F = -F;   
    % end function pdes

% end function shocktube