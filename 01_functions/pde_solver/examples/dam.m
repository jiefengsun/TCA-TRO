function dam
% The breaking dam problem of section 4.5.6 of J. Kevorkian, Partial
% Differential Equations, Analytical Solution Techniques, Wadsworth &
% Brooks/Cole, Pacific Grove, CA, 1990. There is an analytical solution
% of the shallow water equations for quiescent flow that is an expansion
% fan.  When the conservative form of the PDEs is used, the methods compute
% this solution.  There is also an analytical solution with an expansion
% fan and a shock.  Such solutions are worked out in R.J. LeVeque, Finite
% Volume Methods for Hyperbolic Problems, C.U.P., Cambridge, 2002.  When
% the non-conservative form of the PDEs is used, the methods compute this
% solution.  If the initial height h is positive throughout the interval,
% this kind of solution is also computed with the conservative form.
%
% The eigenvalues of a local Jacobian are computed and used in the time
% step function.  They are not distinct when the height h is zero, hence 
% the problem is not strictly hyperbolic then. The height starts out and
% remains zero over part of the interval throughout the integration.
%
% The variables are not the same for the two forms of the PDEs.  In non-
% conservative form the variables are the averaged height h and the
% vertically averaged horizontal speed u.  In conservative form the
% variables are h and u*h.  The initial data specifies h and u = 0 = uh.
% The height variable h is the only one plotted.

global clform

meth = menu('Specify the method','LxF','LxW','NT');
switch meth
    case 1
        method = 'LxF';
    case 2
        method = 'LxW';
    case 3
        method = 'NT';
end
if meth == 3
    clform = 1;
else
    clform = menu('Conservative form of PDEs?','YES','NO');
end
if clform == 1
    form = 3;
    pdefun = @cl;
else
    form = 1;
    pdefun = @pdes;
end

t = 0;
Npoints = 400;
x = linspace(-5,10,Npoints);
dx = x(2) - x(1);
V = zeros(2,Npoints);
V(2,find(x < 0)) = 1;

% Use homogeneous Neumann conditions for quiescent flow.
sol = setup(form,pdefun,t,x,V,method,[],[],{[1,2],[1,2]});

close all
plot(x,V(2,:))
axis([-5 10 -0.1 1.1]) 
title('Initial height of the free surface--dam at x = 0.')
xlabel('t = 0.')
pause

howfar = 1;
for m = 1:4
    if clform == 1
        sol = hpde(sol,howfar,@timestep);
    else
        sol = hpde(sol,howfar,@timestepnc);
    end
    t = sol.t;
    plot(x,sol.u(2,:),'k',x,analytical(t,x),'r')
    legend('Computed','Analytical')
    axis([-5 10 -0.1 1.1]) 
    title(['Height of the free surface at t = ',num2str(t),'.'])
    if clform == 1
        xlabel(['Conservative form of PDEs. ',num2str(Npoints),...
                ' mesh points'])
    else
        xlabel(['Non-conservative form of PDEs. ',num2str(Npoints),...
                ' mesh points'])    
    end
    pause
end

%=========================================================================
% Subfunctions

    function F = cl(V)
        F = zeros(size(V));
        % First component is - (u^2*h + 0.5*h^2).  Recall that V(1,:)
        % provides u*h and V(2,:) provides h.  Must compute carefully 
        % because h starts out and remains 0 over part of the interval.
        F(1,:) = V(1,:).^2;
        ndx = find(V(2,:) > 0);        
        F(1,ndx) = - (F(1,ndx)./V(2,ndx) + 0.5*V(2,ndx).^2);
        F(2,:) = - V(1,:);
    % end function cl
       
    function dt = timestep(dx,t,x,V)
        temp = zeros(size(V(1,:)));
        ndx = find(V(2,:) > 0);
        temp(ndx) = V(1,ndx) ./ V(2,ndx);
        mumax1 = max(abs(temp + sqrt(V(2,:))));
        mumax2 = max(abs(temp - sqrt(V(2,:))));
        dt = 0.9*dx/max(mumax1,mumax2);
    % end function timestep
    
     function F = pdes(t,x,V,V_x)
        F = zeros(size(V));
        F(1,:) = - V(1,:).*V_x(1,:) - V_x(2,:);
        F(2,:) = - V_x(1,:).*V(2,:) - V(1,:).*V_x(2,:);
    % end function pdes
    
    function dt = timestepnc(dx,t,x,V)
        mumax1 = max(abs(V(1,:) + sqrt(V(2,:))));
        mumax2 = max(abs(V(1,:) - sqrt(V(2,:))));
        dt = 0.9*dx/max(mumax1,mumax2);
    % end function timestepnc
        
    function  h = analytical(t,x)
        global clform
        M = length(x);
        h = zeros(1,M);
        if clform == 1
            for m = 1:M
                if x(m) <= -t
                    h(m) = 1;
                elseif x(m) < 2*t
                    h(m) = (1/9)*(2 - x(m)/t)^2;
                end
            end
        else
            hstar = 6 - 4*sqrt(2);  
            ustar = -2 + 2*sqrt(2);
            zetastar = -4 + 3*sqrt(2);
            for m = 1:M
                if x(m) <= -t
                    h(m) = 1;
                elseif x(m) < zetastar*t
                    h(m) = (1/9)*(2 - x(m)/t)^2;
                elseif x(m) < ustar*t
                    h(m) = hstar;
                end
            end
        end
    % end function analytical
                

% end function dam