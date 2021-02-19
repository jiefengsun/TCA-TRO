function sol = hpde3(sol,howfar,timestep)    
% Solve initial-boundary value problem for a hyperbolic system of PDEs 
% of the form u_t = f(u)_x. 
%
% FLUX is the handle of a function to evaluate f(u).  When called with an 
% array U, it is to return an array V with V(:,m) = f(U(:,m)) for each m. 

NeumannL = sol.Neumann{1};
NeumannR = sol.Neumann{2};
have_NL = ~isempty(NeumannL);
have_NR = ~isempty(NeumannR);
flux = sol.pdefun;
bcfun = sol.bcfun;
have_bcfun = ~isempty(bcfun);
variable_step = isa(timestep,'function_handle');

% Local variables
t = sol.t;
x = sol.x;
u = sol.u;

% Get working quantities:
M = length(x);                         % number of mesh points
neqn = size(u,1);                      % number of equations
dx = x(2) - x(1);                      % increment in space

% Allocate working storage
Del = zeros(neqn,M);
D = zeros(neqn,M); 
F = zeros(neqn,M);
Fval = zeros(neqn,M);
fhalf = zeros(neqn,M);
uhalf = zeros(neqn,M);
utemp = zeros(neqn,M);


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

    lambda = dt/dx;   
    half_lambda = lambda/2;
    
    % FIRST HALF STEP
    
    % Approximation of the form considered in Case I:
    % u(:,m) is the cell average for [x(m) - dx/2, x(m) + dx/2]
    % except at the ends where u(:,1) is the cell average for
    % [x(1), x(1) + dx/2] and u(:,M) is the cell average for
    % [x(M) - dx/2, x(M)].  Initially u(:,m) is equal to the
    % given initial value at x(m) for m = 1,...,M.
    
    % Form Del(:,m) = u(:,m+1) - u(:,m).
    Del = diff(u')';  % DIFF does differences by rows; need by columns
    D(:,2:M-1) = MM(Del(:,2:M-1),Del(:,1:M-2));

    % D(:,m) is the scaled slope for u(:,m).  It cannot be defined
    % using the minmod (MM) function for m = 1 and m = M, so linear
    % extrapolation from the interior is used.
    if sol.periodic
        D(:,1) = MM(Del(:,1),u(:,1) - u(:,M-1));
        D(:,M) = D(:,1);
    else
        % Linear extrapolation
        D(:,1) = D(:,2);
        D(:,M) = D(:,M-1);
    end
    
    % Form Del(:,m) = f(:,m+1) - f(:,m).
    F = feval(flux,u);
    Del = diff(F')';  % DIFF does differences by rows; need by columns
    F(:,2:M-1) = MM(Del(:,2:M-1),Del(:,1:M-2));

    % F(:,m) is the scaled flux for f(u(:,m)).  It cannot be defined
    % using the minmod (MM) function for m = 1 and m = M, so linear
    % extrapolation from the interior is used.
    if sol.periodic
        F(:,1) = MM(Del(:,1),F(:,1) - F(:,M-1));
        F(:,M) = F(:,1);
    else
        % Linear extrapolation
        F(:,1) = F(:,2);
        F(:,M) = F(:,M-1);
    end
    
    % uhalf(:,m) is the approximate solution at (x(m),t + dt/4) and
    % fhalf(:,m) is the flux there.
    uhalf(:,1:M) = u(:,1:M) + 0.5*half_lambda*F(:,1:M);
    fhalf(:,1:M) = feval(flux,uhalf);

    % utemp(:,m) approximates the cell average on [x(m),x(m) + dx] at time
    % t + dt/2.
    utemp(:,1:M-1) = 0.5*(u(:,1:M-1) + u(:,2:M)) ...
                     + 0.125*(D(:,1:M-1) - D(:,2:M)) ...
                     + half_lambda*(fhalf(:,2:M) - fhalf(:,1:M-1));
    
    if sol.periodic
        utemp(:,M) = utemp(:,1);
    else
        % Dummy value
        utemp(:,M) = utemp(:,M-1);
    end
    t = t + dt/2;    
    u = utemp;
   
    % Approximation of the form considered in Case II:
    % u(:,m) is the cell average for [x(m), x(m) + dx] for 
    % m = 1,...,M-1.  It is w_{m+1/2}.  
    
    % SECOND HALF STEP
    
    Del = diff(u(:,1:M-1)')';  % Del(:,m) is w_{m+3/2} - w_{m+1/2}
    D(:,2:M-2) = MM(Del(:,2:M-2),Del(:,1:M-3));
    if sol.periodic
        D(:,1) = MM(Del(:,1),u(:,1) - u(:,M-1));
        D(:,M-1) = MM(u(:,1) - u(:,M-1),Del(:,M-2));
    else
        % Linear extrapolation
        D(:,1) = D(:,2);
        D(:,M-1) = D(:,M-2);
    end
    
    % Form Del(:,m) = f(:,m+1) - f(:,m).
    F = feval(flux,u);
    Del = diff(F')';  % DIFF does differences by rows; need by columns
    F(:,2:M-1) = MM(Del(:,2:M-1),Del(:,1:M-2));

    % F(:,m) is the scaled flux for f(u(:,m)).  It cannot be defined
    % using the minmod (MM) function for m = 1 and m = M, so linear
    % extrapolation from the interior is used.
    if sol.periodic
        F(:,1) = MM(Del(:,1),F(:,1) - F(:,M-1));
        F(:,M-1) = MM(F(:,1) - F(:,M-1),Del(:,M-2));
    else
        % Linear extrapolation
        F(:,1) = F(:,2);
        F(:,M) = F(:,M-1);
    end
    uhalf(:,1:M-1) = u(:,1:M-1) + 0.5*half_lambda*F(:,1:M-1);
    fhalf(:,1:M-1) = feval(flux,uhalf(:,1:M-1));

    % utemp(:,m) approximates solution at time t + dt and x(m).
    utemp(:,2:M-1) = 0.5*(u(:,1:M-2) + u(:,2:M-1)) ...
                     + 0.125*(D(:,1:M-2) - D(:,2:M-1)) ...
                     + half_lambda*(fhalf(:,2:M-1)- fhalf(:,1:M-2));
                 
    if sol.periodic  
        utemp(:,1) = 0.5*(u(:,M-1) + u(:,1)) ...
                     + 0.125*(D(:,M-1) - D(:,1)) ...
                     + half_lambda*(fhalf(:,1) - fhalf(:,M-1));
        utemp(:,M) = utemp(:,1);
    else
        uhalfL = (utemp(:,1) - 0.5*D(:,1)) + 0.5*half_lambda*F(:,1);
        fhalfL = feval(flux,uhalfL);
        utemp(:,1) = u(:,1) - 0.25*D(:,1) ...
                     + 2*half_lambda*(fhalf(:,1) - fhalfL);
        uhalfR = (utemp(:,M-1) + 0.5*D(:,M-1)) + 0.5*half_lambda*F(:,M-1);
        fhalfR = feval(flux,uhalfR);
        utemp(:,M) = u(:,M-1) + 0.25*D(:,M-1) ...
                     + 2*half_lambda*(fhalf(:,M-1) - fhalfR);
    end

    
    t = t + dt/2;
    u = utemp;
    
    if ~sol.periodic
        if have_bcfun
            % Apply BCs to cell averages at ends.
            [u(:,1),u(:,M)] = feval(bcfun,t,u(:,1),u(:,M));
        end   
        if have_NL
            u(NeumannL,1) = u(NeumannL,2);
        end
        if have_NR
            u(NeumannR,M) = u(NeumannR,M-1);
        end
    end
    
    sol.nstep = sol.nstep + 1;
    
end

% Load results in structure for return:
sol.t = t;
sol.u = u;

%Subfunction=================================================
    function v = MM(a,b)
        v = 0.5*(sign(a) + sign(b)) .* min(abs(a),abs(b));
    % end function MM    

% end function hpde3

    
    

