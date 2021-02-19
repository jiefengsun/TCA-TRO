function Strikwerda
% Exercises 1.2.3 and 1.2.4 of J.C. Strikwerda, Finite Difference Schemes 
% and Partial Differential Equations, 2nd ed., SIAM, Philadelphia, 2004
% provide the analytical solution of a system of two equations that has a 
% discontinuity because of incompatible initial and boundary data.  There 
% is a Dirichlet boundary condition on one component at each end of the 
% interval.

global A
A = [0 1; 1 0];

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
        source = [0; 0];
        pdefun = {@flux,source};
    case 3
        pdefun = @cl;
end

t = 0;
dx = 0.01;
x = 0:dx:1;
u = ones(2,length(x));
u(1,:) = x;

close all
plot(x,u,'r')
title('Initial data');  
xlabel(['Solution at t = ',num2str(t),'.']);
axis([0,1,-2,2.5])
pause

% A constant step size is appropriate.
timestep = 0.9*dx;

sol = setup(form,pdefun,t,x,u,method,[],@bcfun);

howfar = 0.25;
for m = 1:8
    sol = hpde(sol,howfar,timestep);    
    t = sol.t;
    x = sol.x;
    u = sol.u;
    plot(x,u,'k.-',x,analytical(t,x),'r')
    title(['Computed with ',method,'. Analytical solution in red.'])
    xlabel(['Solution at t = ',num2str(t),'.']);
    axis([0,1,-2,2.5])
    pause   
end

%=========================================================================
% Subfunctions
    
    function [uL,uR] = bcfun(t,uL,uR)
        uL(1) = 0; 
        uR(1) = 0;
    % end function bcfun

    function F = pdes(t,x,u,u_x)
        global A
        F = - A*u_x;
    % end function pdes

    function F = flux(t,x,u)
        global A
        F = - A*u;
    % end function flux
    
    function F = cl(u)
        global A
        F = - A*u;
    % end function flux

    function v = analytical(t,x)
        N = length(x);
        v = zeros(2,N);
        t = min(t,2);
        for m = 1:N
            if t <= 1
                if x(m) < 1 - t
                   v(:,m) = [x(m); 1-t];
                else
                    v(:,m) = [x(m)-1; 2-t];
                end
            else 
                if x(m) < t - 1
                    v(:,m) = [x(m); 3-t];
                else
                    v(:,m) = [x(m)-1; 2-t];
                end
            end
        end
    % end function analytical
            
% end function Strikwerda