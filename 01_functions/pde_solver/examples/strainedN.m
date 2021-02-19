function strainedN
% This program solves numerically the problem that is solved analytically
% in strainedA by the method of strained coordinates.  The PDEs as given
% in the book are u_t + u u_x + h_x = 0 and h_t + u h_x + h u_x = 0.  They
% are solved here in conservative form.  The initial data and the solution 
% are periodic.

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
Npoints = 400;
x = linspace(0,10*pi,Npoints);

% V = [h, u]^T
V = zeros(2,Npoints);
epsilon = 0.1;
V(1,:) = 1 + epsilon*sin(x) + 0.25*epsilon^2*sin(x).^2;  
V(2,:) = epsilon*sin(x);                                 

periodic = true;
sol = setup(3,@cl,t,x,V,method,periodic);

close all
plot(x,V(1,:))
axis([0 10*pi 0 4]) 

hold on
howfar = 2;
for m = 1:4
    sol = hpde(sol,howfar,@timestep);                            
    t = sol.t;
    plot(x,sol.u(1,:)+t/4);
    axis([0 10*pi 0 4]) 
    title(['Breaking wave for \epsilon = ',num2str(epsilon),...
           ' approximated with ',method,'.'])
    xlabel(['Curves give h(x,t) for t = 0:2:8 with vertical offsets of',...
           ' t/4.'])
end
hold off

%=========================================================================
% Subfunctions

    function F = cl(V)
        % h_t = - (u h)_x
        % u_t = - (0.5 u^2 + h)_x 
        F = zeros(size(V));
        F(1,:) = - V(2,:).*V(1,:);
        F(2,:) = - (0.5*V(2,:).^2 + V(1,:));
    % end function cl
    
    function dt = timestep(dx,t,x,V)
        mumax1 = max(abs(V(2,:) + sqrt(V(1,:))));
        mumax2 = max(abs(V(2,:) - sqrt(V(1,:))));
        dt = 0.9*dx/max(mumax1,mumax2);
    % end function timestep               

% end function strainedN