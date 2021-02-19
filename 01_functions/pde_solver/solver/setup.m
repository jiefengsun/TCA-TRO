function sol = setup(form,pdefun,t,x,u,method,periodic,bcfun,Neumann)
% SETUP and HPDE solve initial-boundary value problems for first order
% systems of hyperbolic partial differential equations (PDEs).  The problem
% is defined in SETUP and solved in HPDE.  The PDEs can have three forms:
%
% FORM = 1 Solve u_t = f(t,x,u,u_x). 
%
%          PDEFUN is the handle of a function to evaluate f(t,x,u,u_x).  
%          When called with a scalar T and arrays X,U,U_x, it is to return 
%          an array V with V(:,m) = f(T,X(m),U(:,m),U_x(:,m)) for each m--
%          see NOTE.
%
% FORM = 2 Solve u_t = f(t,x,u)_x + s(t,x,u). 
%
%          The PDEFUN argument has the form {FLUX,SOURCE}. FLUX is the 
%          handle of a function to evaluate f(t,x,u).  When called with a 
%          scalar T and arrays X,U, it is to return an array V with 
%          V(:,m) = f(T,X(m),U(:,m)) for each m--see NOTE. Similarly, SOURCE 
%          is the handle of a function to evaluate s(t,x,u) and it must be 
%          vectorized just like FLUX. If s(t,x,u) is a constant vector, 
%          it can be supplied as SOURCE. 
%
% FORM = 3 Solve u_t = f(u). 
%
%          PDEFUN is the handle of a function to evaluate f(u). When called
%          with an array U, it is to return an array V with V(:,m) = f(U(:,m)) 
%          for each m--see NOTE. 
%
% NOTE     All the input arrays to PDEFUN have the same number of columns, 
%          but this number varies from one call to the next.
%
% The integration starts at time T. The solution is computed on a fixed mesh
% X with initial value U.  For each m, X must satisfy X(m+1) = X(m) + DX for 
% a constant DX > 0 and U(:,m) approximates u(T,X(m)). 
%
% The string METHOD indicates which method is to be used: 
%
%   'LxF'--Lax-Friedrichs method (Shampine's two-step variant).  It is of 
%          order one and dissipative of order two.
%
%   'LxW'--Lax-Wendroff method (Richtmyer's two-step variant).  It is of 
%          order two and dissipative of order four.
%
%   'SLxW'--LxW with a nonlinear filter (Engquist et alia) applied at each
%          time step to reduce the total variation.
%
%   'NT'---Nessyahu-Tadmor method for PDEs of form 3. It is of order two 
%          and total variation diminishing for scalar PDEs.
%
%
% A wide variety of boundary conditions is allowed. Periodic conditions
% cannot be used with other boundary conditions, but homogeneous Neumann 
% and general boundary conditions can be used together.
%
%   Periodic--If the boundary conditions are periodic, set PERIODIC to TRUE 
%             and otherwise, FALSE (or []). 
%
%   General---At each step HPDE computes approximations uLex to u(T,X(1)) 
%             and uRex to u(T,X(end)) by extrapolation from the interior 
%             and passes these vectors to BCFUN along with T. BCFUN is the 
%             handle of a function of the form [uL,uR] = bcfun(t,uLex,uRex)
%             that defines the components of uL and uR so as to impose the
%             boundary conditions at X(1) and X(end), respectively.
%
%   Neumann---If there are homogeneous Neumann boundary conditions at the
%             left end of the interval, i.e. u(T,X(1))_x = 0 for some
%             components I, provide the indices of these components in the
%             vector NeumannL and otherwise set NeumannL = [].  Homogeneous
%             Neumann boundary conditions at the right end of the interval
%             are similarly indicated with NeumannR.  The argument NEUMANN
%             must have the form {NeumannL,NeumannR}. 
%
% To illustrate the coding of boundary conditions, suppose that u(t,x) has 
% 3 components and that the boundary conditions are:  Components 1 and 3 
% are to satisfy homogeneous Neumann conditions at the left end and no
% condition is to be imposed on component 2 there. Component 2 is to have 
% the value sin(t) at the right end and no condition is to be imposed on 
% components 1 and 3 there. These are not periodic boundary conditions, so
% PERIODIC = FALSE.  The homogeneous Neumann conditions are specified by 
% NEUMANN = {[1,3],[]}. Using the same names for both input and output 
% variables allows the remaining conditions to be coded as
%
%   function [uL,uR] = bcfun(t,uL,uR)
%      uR(2) = sin(t);


input_OK = true;
switch form
    case {1,3}
        input_OK = isa(pdefun,'function_handle');
    case 2
        if isa(pdefun,'cell')
            if length(pdefun) == 2
                input_OK = isa(pdefun{1},'function_handle') && ...
                           ( isa(pdefun{2},'function_handle') || ...
                             isnumeric(pdefun{2})  );
            else
                input_OK = false;
            end
        else
            input_OK = false;
        end
end
if ~input_OK
    error('PDEFUN does not agree with FORM.')
end
sol.form = form;
sol.pdefun = pdefun;
        
sol.t = t;
sol.x = x;
sol.u = u;

sol.nstep = 0;

method = upper(method);
if strcmp(method,'NT') && (form ~= 3)
    error('NT can be used only when FORM = 3.')
end
sol.method = method;

switch nargin
    case 6
        sol.periodic = [];
        sol.bcfun = [];        
        sol.Neumann = {[],[]};
    case 7        
        sol.periodic = periodic;
        sol.bcfun = [];        
        sol.Neumann = {[],[]};
    case 8
        sol.periodic = periodic;
        sol.bcfun = bcfun;        
        sol.Neumann = {[],[]};
    case 9
        sol.periodic = periodic;
        sol.bcfun = bcfun;        
        sol.Neumann = Neumann;
end
if isempty(sol.periodic)
    sol.periodic = false;
end
input_OK = true;
if isa(sol.Neumann,'cell')
    input_OK = length(sol.Neumann) == 2;
else
    input_OK = false;
end
if ~input_OK
    error('NEUMANN must have the form {NeumannL,NeumannR}.')
end
