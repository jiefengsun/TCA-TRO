function sol = hpde(sol,howfar,timestep)
% An initial-boundary value problem for a first order system of hyperbolic 
% PDEs and how it is to be solved are defined in SETUP and an initial SOL
% is formed. At each call to HPDE, the input SOL contains fields
%
%   SOL.T--current time t
%   SOL.X--(fixed) spatial mesh x
%   SOL.U--approximation to u(t,X), i.e., for each m, U(:,m) approximates
%          the vector u(t,X(m))
%   SOL.NSTEP--number of steps taken to reach current time
%
% Each call to HPDE advances from SOL.T to SOL.T+HOWFAR.  TIMESTEP has
% two forms:  
%   Positive scalar--solver takes steps of this size. 
%   Handle for a function of the form DT = TIMESTEP(DX,T,X,U)--solver calls 
%        function at each step with current approximation U to the solution 
%        at time T, the (fixed) mesh X, and the constant mesh spacing DX.  
%        The solver takes the step DT returned by TIMESTEP.
% In both cases the last step is reduced as necessary so as to stop exactly
% at SOL.T+HOWFAR. To continue the integration, call HPDE again with this 
% SOL as input. HOWFAR and a scalar TIMESTEP can be changed at each call. 

switch sol.form
    case 1
        sol = hpde1(sol,howfar,timestep);
    case 2
        sol = hpde2(sol,howfar,timestep);
    case 3
        if strcmp(sol.method,'NT')
            sol = hpde3(sol,howfar,timestep);
        else
            sol = hpde2(sol,howfar,timestep);
        end
end

% end function hpde

