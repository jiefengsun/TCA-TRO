function Carleman
% H.G. Kaper and G.K. Leaf study Initial value problems for the Carleman
% equation, Nonlinear Analysis, Theory, Methods & Applications, 4 (1980)
% 343-362.  In particular they consider perturbations about a (positive)
% constant steady state.  They show that convergence to the steady state
% is somewhat faster when the integral of the sum of the two perturbations
% is zero, which is the case here.  The two PDEs arise in form 2.  The
% perturbations are confined to the interval [-pi,pi] and we solve on the
% interval [-infinity,infinity] with infinity = 10.  The steady state
% values are used as boundary conditions at infinity.

meth = menu('Specify the method','LxF','LxW','SLxW');
switch meth
    case 1
        method = 'LxF';
    case 2
        method = 'LxW';
    case 3
        method = 'SLxW';
end

t = 0;
infinity = 10;
Npoints = 400;
x = linspace(-infinity,infinity,Npoints);

V = ones(2,Npoints);
for m = 1:Npoints
    if abs(x(m)) < pi
        V(1,m) = 1 + sin(x(m));
    end
end

close all
plot(x,V(1,:),'b',x,V(2,:)+1.5,'r')
axis([-infinity infinity -0.5 4]) 
title(['Solution at t = ',num2str(t),'.'])
xlabel('Second component (red) plotted with vertical offset of +1.5.')
pause

sol = setup(2,{@flux,@source},t,x,V,method,[],@bcs);

% In the form u_t = A*u_x + S(u), the matrix A is [-1 0; 0 1], so the
% spectral radius is 1.  Accordingly, we use a constant step size.
CFL = 0.9;
dx = x(2) - x(1);
dt = CFL*dx;

% Choose times for output to show the initial behavior and then the
% approach to steady state.
tout = [0 0.05 0.1 0.2 0.5 1 2 4 8 16];
for m = 1:length(tout)-1
    howfar = tout(m+1) - tout(m);
    sol = hpde(sol,howfar,dt);                              
    t = sol.t;
    V = sol.u;
    plot(x,V(1,:),'b',x,V(2,:)+1.5,'r')
    axis([-infinity infinity -0.5 4]) 
    title(['Solution at t = ',num2str(t),' computed with ',method,'.'])
    xlabel('Second component (red) plotted with vertical offset of +1.5.')
    pause
end

%=========================================================================
% Subfunctions

    function F = flux(t,x,V)
        F = [-V(1,:); V(2,:)];
    % end function flux
    
    function S = source(t,x,V)
        temp = V(2,:).^2 - V(1,:).^2;
        S = [temp; -temp];
    % end function source
    
    function [uL,uR] = bcs(t,uL,uR)
        % Steady state values at "infinity".
        uL = [1;1];
        uR = [1;1];
    % end function bcs
              

% end function Carleman