function transmission
% Example 3.3 of Chapter 10 in E.C. Zachmanoglou and D.W. Thoe, Introduction 
% to Partial Differential Equations with Applications, Dover, New York, 1986 
% treats an electrical transmission line.  The PDEs have the form
%   L I_t + E_x + R I = 0
%   C E_t + I_x + G E = 0
% for positive constant parameters L,R,C,G.  The problem is set on the
% whole real line with I(x,0) = f(x), E(x,0) = g(x).  In the special case
% that RC = LG, a distortionless line, an analytical solution is available, 
% so we solve an example of this kind.  We write V = [I,E]^T.  With the 
% support of f(x) and g(x) confined to |x| < 0.5, we solve the problem on
% [-2,2] and impose the undisturbed values of zero as boundary conditions.
%
% The second order LxW does rather better than the first order LxF. There
% is no need for the smoothing of SLxW. Indeed, the filter of SLxW flattens 
% a peak in the plot for t = 0.25 that is present in the analytical solution.

global L R C MF MS
L = 4; R = 4; C = 1; G = 1;
MF = [0, -(1/L); -(1/C), 0];
MS = [-(R/L), 0; 0, -(G/C)];

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
Npoints = 100;
x = linspace(-2,2,Npoints);

V = zeros(2,Npoints);
V(1,:) = f(x);  % I(x,0)
V(2,:) = g(x);  % E(x,0)

close all
plot(x,V,'r')
title('Initial data');  
xlabel(['Solution at t = ',num2str(t),'.']);
axis([-2 2 -0.6 1.1]) 
pause

sol = setup(2,{@flux,@source},t,x,V,method,[],@bcfun);

% The spectral radius of the local Jacobian MF is 1/sqrt(LC),
% so for CFL of 0.9,
dx = x(2) - x(1);
dt = 0.9*dx*sqrt(L*C);  

howfar = 0.25;
for m = 1:5
    sol = hpde(sol,howfar,dt);                              
    t = sol.t;
    V = sol.u;
    plot(x,V,'.-k',x,analytical(t,x),'r')
    axis([-2 2 -0.6 1.1]) 
    title(['Computed with ',method,'. Analytical solution in red.'])
    xlabel(['Solution at t = ',num2str(t),'.']);
    pause
end

%=========================================================================
% Subfunctions

    function v = f(x)
        v = zeros(size(x));
        ndx = find(abs(x) < 0.5);
        v(ndx) = cos(pi*x(ndx));
    % end function f
    
    function v = g(x)
        v = 0.5*f(x);
    % end function g
    
    function F = flux(t,x,V)
        global MF 
        F = MF*V;
    % end function flux

    function S = source(t,x,V)   
        global MS
        S = MS*V;
    % end function source
           
    function [VL,VR] = bcfun(t,VLex,VRex)
        VL = [0;0]; VR = [0;0];
    % end function bcfun
    
    function Vtrue = analytical(t,x)
        global L R C
        fR = f(x - t/sqrt(L*C)); fL = f(x + t/sqrt(L*C));
        gR = g(x - t/sqrt(L*C)); gL = g(x + t/sqrt(L*C));
        Vtrue = zeros(2,length(x));
        Vtrue(1,:) = 0.5*(fR + fL) + 0.5*sqrt(C/L)*(gR - gL);
        Vtrue(2,:) = 0.5*sqrt(L/C)*(fR - fL) + 0.5*(gR + gL);
        Vtrue = Vtrue*exp(-(R/L)*t);      
    % end function analytical
            

% end function transmission