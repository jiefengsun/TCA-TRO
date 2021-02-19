function twostrings
% This is Example 12.6 of A. Stanoyevitch, Introduction to Numerical
% Ordinary and Partial Differential Equations Using MATLAB, Wiley, 
% New York, 2005.  The wave equation with variable sound speed
% u_tt = c^2(t,x,u,u_x) u_xx is used to model a wave propagating through 
% two strings of different density that are fastened together.  The PDE
% is written as a first order system by means of variables v_1 = u, 
% v_2 = u_x, v_3 = u_t.  One of the Dirichlet boundary conditions
% has the form u(0,t) = A(t).  This leads to v_1(0,t) = A(t), but also
% to v_3(0,t) = A'(t).

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
Npoints = 1000;
x = linspace(0,5,Npoints);
% Determine index for mesh point closest to 3:
ndx3 = round((3/5)*Npoints);

% Given u(x,0) = 0, which implies u_x(x,0) = 0, and u_t(x,0) = 0,
% so V is zero at all mesh points when t = 0.
V = zeros(3,Npoints);

sol = setup(1,@pdefun,t,x,V,method,[],@bcfun);

% The spectral radius of the local Jacobian is |c(x)|. The maximum over x
% is 2, so for CFL of 0.9,
dx = x(2) - x(1);
dt = 0.9*dx/2;  

close all
plot(x,sol.u(1,:),'b',x(ndx3),sol.u(1,ndx3),'ro')
title(['At t = ',num2str(sol.t),'.'])
axis([0 5 -1.1 1.1])
xlabel('Two strings joined at red circle.')
pause

howfar = 0.5;
for count = 2:11
    sol = hpde(sol,howfar,dt); 
    plot(x,sol.u(1,:),'b',x(ndx3),sol.u(1,ndx3),'ro')
    title(['At t = ',num2str(sol.t),'.'])
    axis([0 5 -1.1 1.1])
    xlabel('Two strings of different density joined at red circle.')
    pause
end

%=========================================================================
% Subfunctions

    function F = pdefun(t,x,V,V_x)
        F = zeros(size(V));
        F(1,:) = V(3,:);
        F(2,:) = V_x(3,:);
        % c(x) = 2 if 0 <= x <= 3, 1 otherwise. csq is c^2(x).
        csq = ones(size(x));
        csq(find(x <= 3)) = 4;
        F(3,:) = csq.*V_x(2,:);
    % end function pdefun
    
    function [VL,VR] = bcfun(t,VL,VR)
        if t <= pi/5
            VL(1) = sin(5*t);
            VL(3) = 5*cos(5*t);
        else
            VL(1) = 0;
            VL(3) = 0;
        end
        VR(1) = 0;
        VR(3) = 0;
    % end function bcfun

% end function twostrings