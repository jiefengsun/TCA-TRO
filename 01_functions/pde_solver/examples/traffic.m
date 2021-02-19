function traffic
% Problem B.9 from K.E. Gustafson, Partial Differential Equations and 
% Hilbert Space Methods, 2nd ed., Wiley, New York, 1987. This problem is
% an instance of a model of traffic flow discussed on p. 266 ff.  He reports 
% that the leapfrog method and two codes based on the Lax-Wendroff method 
% "... failed to produce adequate backward movement of the traffic bulge."  
% The figure on p. 393 was computed using the method of characteristics.  
% Gustafson presents the PDE in non-conservative form, but does say that 
% Lax-Wendroff schemes often treat problems in conservative form.  He does 
% not discuss boundary conditions nor provide any details about the mesh. 
% Here we use a homogeneous Neumann condition at x = 1 and no boundary 
% condition at x = 0. This interval is big enough that the bulge does not 
% reach a boundary by t = 0.009. A rather fine mesh seems to be necessary 
% to resolve adequately the traffic density at this time, which appears to 
% be the time at which a shock forms.  

meth = menu('Specify the method','LxF','LxW','SLxW'); % provide a user interface
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
x = linspace(0,1,Npoints);
rho0 = 200 + 100*exp(-0.5*exp(1)*((x - 0.5)/0.04).^2);

% Homogeneous Neumann boundary condition at the right end:
sol = setup(1,@pdefun,t,x,rho0,method,[],[],{[],[1]});
sol = hpde(sol,0.009,@timestep);

fprintf('The integration required %i steps.\n',sol.nstep)

close all
plot(x,rho0,'b',x,sol.u,'r')
axis([0 1 150 350])
legend('At t = 0.','At t = 0.009.')
title(['Solution computed with ',method,' using ',...
       num2str(Npoints),' equally spaced mesh points.'])
   
%=========================================================================
% Subfunctions
    
    function F = pdefun(t,x,rho,rho_x)
        F = -c(rho).*rho_x;
    % end function pdefun
    
    function v = c(rho)
        v = 76.184 - 17.2*log(rho);
    % end function c
    
    function dt = timestep(dx,t,x,rho)
        dt = 0.9*dx/max(abs(c(rho)));
    % end function timestep
        
    
% end function traffic
