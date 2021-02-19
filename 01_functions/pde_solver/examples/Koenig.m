function Koenig
% Problem from D.M. Koenig, Invariant imbedding: new design method in unit
% operations, Chem. Engrg., 74 (1967) 181-184. There are a number of steady
% states, but the one computed here incorporates the equilibrium data.

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
Npoints = 200;
x = linspace(0,10,Npoints);
u = zeros(1,Npoints);

sol = setup(1,@pdes,t,x,u,method,[],@bcs);
tout = [0 0.05 0.1 0.3 0.8 1.5];

close all
hold on
for i = 1:length(tout)-1
    howfar = tout(i+1) - tout(i);
    sol = hpde(sol,howfar,@timestep);
    plot(sol.x,sol.u)
    axis([0 10 0 5])
end
plot(x,steady(x),'r')
title(['Solutions for t = 0.05, 0.1, 0.3, 0.8, 1.5 computed with ',...
        method,'.'])
xlabel('One analytical steady state is solid line in red.')
hold off

%=========================================================================
% Subfunctions
    
    function F = pdes(t,x,u,u_x)
        F = -A(t,x,u).*u_x + S(t,x,u);
    % end function pdes
    
    function a = A(t,x,u)
        a = 2*(x - 2*u - 0.01*u.^2).*(1 - x);
    % end function A
    
    function s = S(t,x,u)
        s = 2*(x - 2*u - 0.01*u.^2).*(1 - u);
    % end function S
    
    function [uL,uR] = bcs(t,uL,uR)
        uL = 0;
    % end function bcs

    function v = steady(x)
        % There are several steady states.  One is 
        % found by solving 0.01*u^2 + 2*u - x = 0.
        v = (-2 + sqrt(4 + 4*0.01*x))/(2*0.01);
    % end function steady
    
    function dt = timestep(dx,t,x,u)
        dt = 0.9*dx/max(abs(A(t,x,u)));
    % end function timestep
        
    
% end function Koenig
