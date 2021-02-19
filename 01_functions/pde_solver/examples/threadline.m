function threadline
% Problem corresponding to Fig. 1 of R.D. Swope and W.F. Ames, Vibrations 
% of a moving threadlne, J. Franklin Inst., 275 (1963) 36-55.  They write
% the PDEs as a pair of first order equations in a vector W = [y,u].  Only 
% y is of interest.  Swope & Ames discuss various values of a parameter 
% alpha, of which +7/4 is an interesting choice because of a kink that
% develops in the solution.

global A

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

alpha = +7/4;
beta = 0.25*alpha^2 - 1;
A = [0 1; -beta alpha];

t = 0;
Npoints = 200;
x = linspace(0,1,Npoints);

W = zeros(2,Npoints);
W(1,:) = 0.1*sin(pi*x);

close all
subplot(3,2,1), plot(x,W(1,:))
axis([0 1 -.07 .15]) 
xlabel('t = 0.')

sol = setup(3,@cl,t,x,W,method,[],@bcs);

% A constant step size is appropriate.
mu = max(abs(alpha/2 + 1),abs(alpha/2 - 1));
dx = x(2) - x(1);
dt = 0.9*dx/mu;

howfar = 0.5;
for i = 1:4
    sol = hpde(sol,howfar,dt);
    t = sol.t;
    W = sol.u;
    subplot(3,2,i+1),plot(x,W(1,:))
    axis([0 1 -.07 .15])
    xlabel(['t = ',num2str(i/2),'.'])
end

%=========================================================================
% Subfunctions

    function F = cl(u)
        global A
        F = -A*u;
    % end function cl

    function [uL,uR] = bcs(t,uL,uR)
        uL(1) = 0;
        uR(1) = 0.1*sin(t*pi/2);
    % end function bcs

% end function threadline