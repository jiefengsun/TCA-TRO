function sol = hpde2(sol,howfar,timestep)
% Solve initial-boundary value problem for a hyperbolic system of PDEs 
% of the form u_t = f(t,x,u)_x + s(t,x,u). 
%
% FLUX is the handle of a function to evaluate f(t,x,u).  When called 
% with a scalar T and arrays X,U, it is to return an array V with 
% V(:,m) = f(T,X(m),U(:,m)) for each m. Similarly, SOURCE is the handle 
% of a function to evaluate s(t,x,u) and it must be vectorized just like 
% FLUX. If s(t,x,u) is a constant vector, the vector can be supplied as 
% SOURCE.
%
% This program is also used to solve conservation laws u_t = f(u)_x  
% (FORM = 3) with the LxF and LxW methods as a special case of FORM = 2. 

global pass_flux

method = upper(sol.method);
if strcmp('LXF',method)
    LxF = true;
    smooth = false;
elseif strcmp('LXW',method)
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

if sol.form == 2
    flux = sol.pdefun{1};
    source = sol.pdefun{2};
else
    pass_flux = sol.pdefun;
    flux = @fluxCL;
    source = zeros(npdes,1);
end
constant_S = isnumeric(source);
if constant_S
    Sval = repmat(source(:),1,M-1);
else
    Sval = zeros(npdes,M-1);
end

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
    Fval = feval(flux,t,x,u);
    umid  = (u(:,2:M) + u(:,1:M-1))/2;
    if ~constant_S
        Sval = feval(source,t,xmid,umid);
    end
    vhalf(:,1:M-1) = umid + ...
                        0.5*dt*( Sval + (Fval(:,2:M) - Fval(:,1:M-1))/dx );
    t = t + 0.5*dt;

    % Half step with LxF or full step with leapfrog for LxW:
    Fval = feval(flux,t,xmid,vhalf(:,1:M-1));
    umid(:,2:M-1)  = (vhalf(:,2:M-1) + vhalf(:,1:M-2))/2;
    if ~constant_S
        Sval(:,2:M-1) = feval(source,t,x(2:M-1),umid(:,2:M-1));
    end

    if sol.periodic
        utemp =  0.5*(vhalf(:,1) + vhalf(:,M-1));
        if ~constant_S
            Stemp = feval(source,t,x(1),utemp);
        else
            Stemp = source;
        end
        if LxF
            u(:,1) = utemp + 0.5*dt*( Stemp + (Fval(:,1) - Fval(:,M-1))/dx );
        else
            u(:,1) = u(:,1) + dt*( Stemp + (Fval(:,1) - Fval(:,M-1))/dx );
        end 
        u(:,M) = u(:,1);
    end
    
    if LxF
        u(:,2:M-1) = umid(:,2:M-1) + ...
             0.5*dt*( Sval(:,2:M-1) + (Fval(:,2:M-1) - Fval(:,1:M-2))/dx );
    else
        u(:,2:M-1) = u(:,2:M-1) + ...
                 dt*( Sval(:,2:M-1) + (Fval(:,2:M-1) - Fval(:,1:M-2))/dx );
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

%Subfunction=================================================
    function v = fluxCL(t,x,u)
        global pass_flux
        v = feval(pass_flux,u);
    % end function fluxCL   
    
% end function hpde2

    
    

