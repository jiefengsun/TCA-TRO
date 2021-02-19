function waveeqn
% Example 7.1 of Chapter 8 in E.C. Zachmanoglou and D.W. Thoe, Introduction 
% to Partial Differential Equations with Applications, Dover, New York, 1986 
% treats a semi-infinite string with fixed end. Fig. 7.4 shows how an initial
% profile splits into waves travelling both left and right.  When the wave
% moving left reaches x = 0, it is reflected and starts moving right. Here 
% the problem is solved on [0,10] with an extrapolation boundary condition 
% at the right end.  The wave equation is written as a first order system 
% using the variables of B. Gustafsson, H.-O. Kreiss, and J. Oliger, Time 
% Dependent Problems and Difference Methods, Wiley, New York, 1995, found 
% on p. 364. The PDEs are written as a conservation law so that all the 
% methods can be applied.

t = 0;
infinity = 10;
Npoints = 400;
x = linspace(0,infinity,Npoints);

V = zeros(2,Npoints);
% The desired solution u corresponds to V(1,:). The other variable v is of
% no interest except that v(x,0) is the integral from 0 to x of u_t(x,0). 
% With the given value u_t(x,0) = 0, we have V(2,:) = 0.
for m = 1:Npoints
    if abs(x(m) - 3) < 0.5
        V(1,m) = cos(pi*(x(m)-3));
    end
end

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

close all
plot(x,V(1,:))
axis([0 infinity -0.6 1.1]) 
title('Initial data will split into waves moving in both directions.')
xlabel('t = 0.')
pause

% A constant step size is appropriate.
dx = x(2) - x(1);
timestep = 0.9*dx;

sol = setup(3,@cl,t,x,V,method,[],@bcfun);

howfar = 1;
for m = 1:6
    sol = hpde(sol,howfar,timestep);                              
    t = sol.t;
    V = sol.u;
    plot(x,V(1,:))
    axis([0 infinity -0.6 1.1]) 
    xlabel(['t = ',num2str(t),'.'])
    if m <= 3
        title('Waves moving both left and right.')
    else
        title('Wave moving left reflects at x = 0 and starts moving right.')
    end
    pause(1)
end

%=========================================================================
% Subfunctions

    function F = cl(V)
        F = [V(2,:); V(1,:)];
    % end function cl
    
    function [VL,VR] = bcfun(t,VL,VR)
        % u(0,t) = 0, extrapolation BC at infinity
        VL(1) = 0;
    % end function bcfun

% end function waveeqn