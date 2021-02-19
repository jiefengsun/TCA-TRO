function ionx
% Equilibrium ion exchange problem from S. Goldstein, On the mathematics
% of exchange processes in fixed columns. II, Proc. Roy. Soc., 219 (1953)
% 171-185.  This is a special case that Goldstein solves analytically.
% The PDE is a conservation law u_t + g(u)_x = 0, though the analysis is
% done on the non-conservative form u_t + g'(u)u_x = 0. When the parameter
% r > 1, the solution has two kinds of behavior.  A discontinuity and 
% discontinuous derivatives present numerical difficulties.  The second
% order methods deal better with one of the discontinuous derivatives. The
% nonlinear filter of 'SLxW' is quite helpful.

global r X Npoints

% There are two parameters
r = 2;
X = 1;

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
Npoints = 600;
x = linspace(0,4,Npoints);
u0 = zeros(size(x));
u0(find(x <= X)) = 1;
u0(1) = 0;

close all
plot(x,u0)
title('Initial data')
axis([-0.1 4 -0.1 1.1])
pause

sol = setup(3,@cl,t,x,u0,method,[],@bcfun);

% The behavior of the solution changes at t = X/(r-1).
howfar = 0.2*(2*X/(r-1));
for i = 1:5
    sol = hpde(sol,howfar,@timestep);
    t = sol.t;
    u = sol.u;
    v = analytical(t,x);
    plot(x,u,'b',x,v,'r')
    axis([-0.1 4 -0.1 1.1])
    title(['Computed with ',method,' using ',...
           num2str(Npoints),' equally spaced mesh points.'])
    xlabel(['Solution at t = ',num2str(sol.t),'.']);
    legend('Computed','Analytical')
    pause
end

%=========================================================================
% Subfunctions
    
    function F = cl(u)
        global r
        F = - u ./ (r + (1-r)*u);
    % end function cl
    
    function dt = timestep(dx,t,x,u)
        global r
        gp = r ./ (r + (1-r)*u).^2;
        dt = 0.9*dx/max(abs(gp));
    % end function timestep
    
    function [uL,uR] = bcfun(t,uL,uR)
        uL(1) = 0;
        uR(1) = 0;
    % end function bcfun            
    
    function v = analytical(t,x)
        global r X Npoints
        v = zeros(1,Npoints);
        if t <= X/(r-1)
            for i = 1:Npoints
                if x(i) <= t/r
                    v(i) = 0;
                elseif x(i) <= t*r
                    v(i) = (r - sqrt(r*t/x(i)))/(r-1);
                elseif x(i) <= X + t
                    v(i) = 1;
                else
                    v(i) = 0;
                end
            end
        else
            xstar = (sqrt(t) + sqrt((r-1)*X))^2 / r;
            for i = 1:Npoints
                if x(i) <= t/r
                    v(i) = 0;
                elseif x(i) <= xstar
                    v(i) = (r - sqrt(r*t/x(i)))/(r-1);
                else
                    v(i) = 0;
                end
            end
        end
            
    % end function analytical
                    
    
% end function ionx
