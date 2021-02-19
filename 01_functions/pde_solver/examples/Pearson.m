function Pearson
% C.E. Pearson, Dual time scales in a wave problem governed by coupled
% nonlinear equations, SIAM Review, 23 (1981) 425-433, approximates a
% periodic solution of the shallow water equations by perturbation and
% direct numerical solution and plots it in Fig. 1.  The PDEs as given
% as u_t + u u_x + n_x = 0 and n_t + [u*(1 + n)]_x = 0.  They are solved 
% here in conservative form. Pearson studies two problems, one with periodic 
% data on [0,1] and one with an isolated disturbance on (-infinity,infinity).  
% Both are solved in this example. The periodic problem is solved easily  
% by all the methods with a modest number of mesh points.  The isolated 
% disturbance is confined to |x| < 1/2.  The discontinous disturbance splits 
% and moves in both directions, reaching about -15 and +15 in the time 
% specified.  The LxW method provides an acceptable solution, but it has an 
% oscillation that is not physical.  Filtering improves this. The first order 
% LxF is qualitatively correct, but the solution is rather damped.  The second 
% order NT provides the best solution.

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

% Periodic disturbance
t = 0;
Npoints = 100;
x = linspace(0,1,Npoints);
% V = [u,n]^T
V = zeros(2,Npoints);
epsilon = 0.01;
V(2,:) = epsilon*cos(2*pi*x);   

periodic = true;
sol1 = setup(3,@cl,t,x,V,method,periodic);
T = 15;
sol1 = hpde(sol1,T,@timestep);

close all
plot(sol1.x,sol1.u(2,:))
axis([0 1 -0.01 +0.01])
title(['Surface elevation \eta(x,',num2str(T),') computed with ',method,'.'])
xlabel(['Periodic disturbance with \epsilon = ',num2str(epsilon),'.'])

% Isolated initial disturbance
t = 0;
infinity = 20;
Npoints = 2000;
x = linspace(-infinity,infinity,Npoints);
V = zeros(2,Npoints);
epsilon = 0.01;
ndx = find(abs(x) < 0.5);
V(2,ndx) = epsilon*(1 + cos(2*pi*x(ndx)));

sol2 = setup(3,@cl,t,x,V,method,[],@bcs);
T = 15;
sol2 = hpde(sol2,T,@timestep);

figure
plot(sol2.x,sol2.u(2,:))
axis([-infinity infinity -0.005 +0.015])
title(['Surface elevation \eta(x,',num2str(T),') computed with ',method,'.'])
xlabel(['Isolated initial disturbance with \epsilon = ',num2str(epsilon),'.'])

figure
plot(sol2.x,sol2.u(2,:))
axis([14 16 -0.005 +0.015])
title(['Surface elevation \eta(x,',num2str(T),') computed with ',method,'.'])
xlabel(['Isolated initial disturbance with \epsilon = ',num2str(epsilon),'.'])

%=========================================================================
% Subfunctions

    function F = cl(V)
        % V = [u,n]^T
        % u_t = - (0.5 u^2 + n)_x
        % n_t = - [u*(1 + n)]_x    
        F = zeros(size(V));
        F(1,:) = - (0.5*V(1,:).^2 + V(2,:));        
        F(2,:) = - V(1,:).*(1 + V(2,:));
    % end function cl
    
    function dt = timestep(dx,t,x,V)
        mumax1 = max(abs(V(1,:) + sqrt(1+V(2,:))));
        mumax2 = max(abs(V(1,:) - sqrt(1+V(2,:))));
        dt = 0.9*dx/max(mumax1,mumax2);
    % end function timestep    
    
    function [uL,uR] = bcs(t,uLex,uRex)
        uL = [0;0]; 
        uR = [0;0];
    % end function bcs
            
% end function Pearson