function sol = hpde1(sol,howfar,timestep)
% Solve initial-boundary value problem for a hyperbolic system of PDEs 
% of the form u_t = f(t,x,u,u_x). 
%
% PDEFUN is the handle of a function to evaluate f(t,x,u,u_x).  When 
% called with a scalar T and arrays X,U,U_x, it is to return an array V 
% with V(:,m) = f(T,X(m),U(:,m),U_x(:,m)) for each m.

if strcmp('LXF',sol.method)
    LxF = true;
    smooth = false;
elseif strcmp('LXW',sol.method)
    LxF = false;
    smooth = false;
else
    LxF = false;
    smooth = true;
end
NeumannL = sol.Neumann{1};
NeumannR = sol.Neumann{2};
have_NL = ~isempty(NeumannL);
have_NR = ~isempty(NeumannR);
bcfun = sol.bcfun;
if isempty(bcfun)
    bcfun = @dummyBC;
end
pdefun = sol.pdefun;
variable_step = isa(timestep,'function_handle');

% Local variables
t = sol.t;
x = sol.x;
u = sol.u;

M = length(x);                         % number of mesh points
dx = x(2) - x(1);                      % increment in space
npdes = size(u,1);                     % number of PDEs
vhalf = zeros(npdes,M);                % working storage

xmid = (x(2:M) + x(1:M-1))/2;
tout = t + howfar;
done = false;
while ~done
    
    if variable_step
        dt = feval(timestep,dx,t,x,u);
    else
        dt = timestep;
    end
    if t+dt >= tout
        dt = tout - t;
        done = true;
    end 

    % Half step with LxF.
    umid  = (u(:,2:M) + u(:,1:M-1))/2;
    uxmid = (u(:,2:M) - u(:,1:M-1))/dx;
    Fval = feval(pdefun,t,xmid,umid,uxmid);
    vhalf(:,1:M-1) = umid + 0.5*dt*Fval;
    t = t + 0.5*dt;

    % Half step with LxF or full step with leapfrog for LxW:
    if sol.periodic
        pxmid = x(1);
        pumid =  (vhalf(:,1) + vhalf(:,M-1))/2;
        puxmid = (vhalf(:,1) - vhalf(:,M-1))/dx;
        pFval = feval(pdefun,t,pxmid,pumid,puxmid);
        if LxF
            u(:,1) = pumid + 0.5*dt*pFval;
        else
            u(:,1) = u(:,1) + dt*pFval;
        end 
        u(:,M) = u(:,1);
    end
    
    umid  = (vhalf(:,2:M-1) + vhalf(:,1:M-2))/2;
    uxmid = (vhalf(:,2:M-1) - vhalf(:,1:M-2))/dx;
    Fval = feval(pdefun,t,x(2:M-1),umid,uxmid);
    if LxF
        u(:,2:M-1) = umid + 0.5*dt*Fval;
    else
        u(:,2:M-1) = u(:,2:M-1) + dt*Fval;
    end
    t = t + 0.5*dt;
    
    if ~sol.periodic
        % Extrapolate to form approximations uLex and uRex:
        uLex = 2*u(:,2)   - u(:,3);
        uRex = 2*u(:,M-1) - u(:,M-2);
        
        [u(:,1),u(:,M)] = feval(bcfun,t,uLex,uRex);
        
        if have_NL
            u(NeumannL,1) = u(NeumannL,2);
        end
        if have_NR
            u(NeumannR,M) = u(NeumannR,M-1);
        end
    end
    
    if smooth
        for eqn = 1:npdes
            u(eqn,1:M) = nlfilter2p1(u(eqn,1:M));
        end
    end  
    
    sol.nstep = sol.nstep + 1;

end

% Load results in structure for return:
sol.t = t;
sol.u = u;

% end function hpde1

