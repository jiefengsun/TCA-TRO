function Schiesser
% W.E. Schiesser, Computational Mathematics in Engineering and Applied
% Science, CRC Press, Boca Raton, FL, 1994, P. 195 ff develops a model of
% heat transfer in section 3.1. He works out semi-analytical solutions for 
% four sets of initial-boundary data.  An incompatibility of initial and
% boundary data for some sets leads to a discontinuity that propagates from
% left to right.  Profiles of the solution are plotted for t = 0.1:0.1:0.5
% to show this.  A surface plot using values computed for t = 1:1:10 show
% the overall behavior of the solution.  Schiesser uses only 21 equally
% spaced mesh points, but 100 are used here to obtain smooth profiles. The
% values he reports for u_1(1,t) at t = 1:1:10 are used to measure the
% error of the numerical solution.

global c1 c2 v Scase
c1 = 0.56; c2 = 0.1; v = 2.031;

Scase = menu('Specify the data set','1','2','3','4');           
 
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
Nx = 100;
x = linspace(0,1,Nx);
u = zeros(2,Nx);

sol = setup(2,{@flux,@source},t,x,u,method,[],@bcs);

% In the form u_t = (A*u)_x + S(u), the matrix A is [-v 0; 0 0], so the
% spectral radius is v (assumed positive).  Accordingly, we use a constant 
% step size with a CFL number of 0.9.
dx = x(2) - x(1);
dt = 0.9*dx/v;

tout = [0:0.1:0.5, 1:1:10];
Nt = length(tout);
% u1 holds output values for u_1(x,t).  It is zero for t = 0.
u1 = zeros(Nt,Nx);  
for m = 2:Nt
    howfar = tout(m) - tout(m-1);
    sol = hpde(sol,howfar,dt);
    u1(m,:) = sol.u(1,:);
end

close all
plot(x,u1(2:6,:))
title(['u_1(x,t) for t = 0.1:0.1:0.5 computed with ',method,'.'])
switch Scase
    case 1
        xlabel(['Data set 1: u_1(0,t) = 1 - exp(-',num2str(c2),'*t)'])
    case 2
        xlabel(['Data set 2: u_1(0,t) = exp(-',num2str(c2),'*t)'])
    case 3
        xlabel('Data set 3: u_1(0,t) = 1')
    case 4
        xlabel('Data set 4: u_1(0,t) = min(t,1)')
end

figure
surf(x,1:10,u1(7:end,:))
xlabel(['Distance x'])
ylabel('Time t')
title(['u_1(x,t) for data set ',num2str(Scase),' computed with ',...
       method,' using ',num2str(Nx),' mesh points.'])

% Reference values at x = 1 and t = 1:1:10 for the four data sets
% taken from Schiesser's book.
u1true = zeros(4,10);
u1true(1,:) = [0.0379, 0.1085, 0.1741, 0.2350, 0.2916, ...
               0.3442, 0.3931, 0.4384, 0.4804, 0.5194];
u1true(2,:) = [0.7316, 0.6802, 0.6322, 0.5874, 0.5456, ...    
               0.5066, 0.4702, 0.4364, 0.4048, 0.3754];
u1true(3,:) = [0.7695, 0.7887, 0.8063, 0.8224, 0.8372, ...
               0.8508, 0.8633, 0.8747, 0.8852, 0.8948];
u1true(4,:) = [0.3880, 0.7792, 0.7975, 0.8144, 0.8299, ...
               0.8441, 0.8571, 0.8690, 0.8800, 0.8900];
max_err = norm(u1true(Scase,:) - u1(7:end,end)');
fprintf('Maximum error in u_1(1,t) for t = 1:10 was %e.\n',max_err) 

%=========================================================================
% Subfunctions

    function F = flux(t,x,u)
        global c1 c2 v Scase
        F = zeros(size(u));
        F(1,:) = -v*u(1,:);
    % end function flux
    
    function S = source(t,x,u)
        global c1 c2 v Scase
        temp = u(2,:) - u(1,:);
        S = [c1*temp; -c2*temp];
    % end function source
    
    function [uL,uR] = bcs(t,uL,uR)
        global c1 c2 v Scase
        switch Scase
            case 1
                uL(1) = 1 - exp(-c2*t);
            case 2
                uL(1) = exp(-c2*t);
            case 3
                uL(1) = 1;
            case 4
                uL(1) = min(t,1);
        end
    % end function bcs
            
    
% end function Schiesser