function Davis
% Burgers' equation as solved in S.F. Davis, A simplified TVD finite 
% difference scheme via artificial viscosity, SIAM J. Sci. Stat. Comput. 
% 6 (1987) 1-18. The problem is solved with periodic boundary conditions. 
% The initial data is a square wave, but the interval is not stated, nor 
% are the details of the square wave. Davis observes that the solution 
% depends strongly on the initial data. Here the square wave is given both 
% negative and positive values so that the upwind direction changes in the 
% span of the interval. Though not the same initial-boundary value problem, 
% the results of this program can be compared qualitatively to Figures 3a-3c 
% of Davis' paper.

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
        source = 0;
        pdefun = {@flux,source};
    case 3
        pdefun = @cl;
end

t = 0;
Npoints = 100;
x = linspace(-2,2,Npoints);
dx = x(2) - x(1);
u = zeros(1,Npoints);
for i = 1:Npoints
    if abs(x(i)) <= 1/2
        u(i) = 1/2;
    else
        u(i) = -1/4;
    end
end

close all
plot(x,u,'r.-')
title('Initial data'); 
xlabel(['Solution at t = ',num2str(t),'.']);
axis([-2 2 -1 1])
pause

periodic = true;
sol = setup(form,pdefun,t,x,u,method,periodic);
  
howfar = 1;
for m = 1:4
    sol = hpde(sol,howfar,@timestep);
    t = sol.t;
    x = sol.x;
    u = sol.u;
    plot(x,u,'r.-')
    title(['Computed with ',method,'.'])
    xlabel(['Solution at t = ',num2str(t),'.']);
    axis([-2 2 -1 1])
    pause    
end

%=========================================================================
% Subfunctions

    function F = pdes(t,x,u,u_x)
        F = - u.*u_x;
    % end function pdes
    
    function v = flux(t,x,u)
        v = -0.5*u.^2;
    % end function flux
    
    function v = cl(u)
        v = -0.5*u.^2;
    % end function cl
    
    function dt = timestep(dx,t,x,u)
        dt = 0.9*dx/max(abs(u));
    % end function timestep
           
% end function Davis