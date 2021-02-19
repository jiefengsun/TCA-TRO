function weakinstab
% Weak nonlinear instability problem. A. Nayfeh, Perturbation Methods,
% Wiley, New York, 1973, solves u_{tt} - u_{xx} - u = u^3 with initial
% conditions u(x,0) = epsilon*cos(K*x), u_t(x,0) = 0 for small epsilon.
% Nayfeh uses the method of renormalization to get a uniformly valid
% asymptotic approximation to the solution when K > 1. (Modes with K < 1
% are unstable.) Here we solve the problem for epsilon = 0.1 and K = 3.  
% G.B. Whitham, Linear and Nonlinear Waves, Wiley, New York, 1974, uses 
% this PDE without its nonlinear term to illustrate two ways of writing a 
% second order equation as a system of first order equations. One is to use 
% variables V(1) = u, V(2) = u_t, V(3) = u_x. Only V(1) is of interest. The 
% problem is solved with periodic boundary conditions and the numerical 
% solution is compared to the asymptotic approximation.  For this problem 
% smoothing is not helpful and when it is used in an integration to times
% that are O(1/epsilon), the solver begins tracking an unstable mode.

global A epsilon sigma K

meth = menu('Specify the method','LxF','LxW','SLxW');
switch meth
    case 1
        method = 'LxF';
    case 2
        method = 'LxW';
    case 3
        method = 'SLxW';
end
form = menu('Specify the form of the PDEs','1','2');
if form == 1
    pdefun = @pdes;
else
    pdefun = {@flux,@source};
end

A = [0 0 0; 0 0 1; 0 1 0];
epsilon = 0.1;
K = 3;
sigma = sqrt(K^2 - 1)*(1 - 9*epsilon^2/(32*(K^2 - 1)));

t = 0;
Npoints = 100;
x = linspace(-pi,pi,Npoints);
dx = x(2) - x(1);
V = zeros(3,Npoints);
V(1,:) = epsilon*cos(K*x);        % u(0,x)
V(2,:) = 0;                       % u_t(0,x)
V(3,:) = -K*epsilon*sin(K*x);     % u_x(0,x)

% The spectral radius of A is 1, so for CFL of 0.9,
dt = 0.9*dx;  

periodic = true;
sol = setup(form,pdefun,t,x,V,method,periodic);

close all
plot(x,V(1,:))
axis([-pi pi -2*epsilon +2*epsilon]) 
title(['With \epsilon = ',num2str(epsilon),' and K = ',num2str(K),'.'])
xlabel('t = 0.')
pause

howfar = 1;
for m = 1:6 
    sol = hpde(sol,howfar,dt); 
    t = sol.t;
    V = sol.u;
    plot(x,V(1,:),'r.',x,asymptotic(t,x),'k')
    title(['With \epsilon = ',num2str(epsilon),' and K = ',num2str(K),'.'])
    legend(method,'Asymptotic approximation')
    axis([-pi pi -2*epsilon +2*epsilon]) 
    xlabel(['t = ',num2str(t),'.'])
    pause
end

%=========================================================================
% Subfunctions

    function F = pdes(t,x,V,V_x)
        global A
        F = A*V_x + source(t,x,V);
    % end function pdes

    function F = flux(t,x,V)
        global A
        F = A*V;
    % end function flux

    function s = source(t,x,V)
        s = zeros(size(V));
        s(1,:) = V(2,:);
        s(2,:) = V(1,:).*(1 + V(1,:).^2);
    % end function source
    
    function u = asymptotic(t,x)
        global epsilon sigma K
        u = epsilon*cos(sigma*t)*cos(K*x);
    % end function asymptotic

% end function weakinstab