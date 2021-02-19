function strainedA
% E.J. Hinch, Perturbation Methods, C.U.P., Cambridge, 1991, studies the
% breaking of shallow water waves by the method of strained coordinates
% in Section 6.3.  This program computes his Fig. 6.3.  He shows that the
% wave breaks at time t = (2/3)(1/epsilon) + O(epsilon^2).  Here epsilon is
% 0.1, so it breaks at about t = 6.67.  The figure shows the wave at times
% t = 6 and 8.  Note that the curve is multivalued at t = 8. For comparison
% with the numerical results, it is more convenient to plot on [0,10*pi]
% than on the interval [0,30] of Fig. 6.3.

epsilon = 0.1;
beta = linspace(-10,40,200);
close all
hold on
for t = 0:2:8
    alpha = getalpha(t,epsilon,beta);
    x = 0.5*(alpha+beta) + epsilon*( 0.3750*(alpha-beta).*sin(beta) ...
        + 0.1250*(cos(beta) - cos(alpha)));
    h = 1 + epsilon*sin(beta) + 0.25*epsilon^2*sin(beta).^2;
    plot(x,h+0.25*t)
    axis([0 10*pi 0 4])
    title(['Breaking wave for \epsilon = ',num2str(epsilon),...
           ' approximated by method of strained coordinates.'])
    xlabel(['Curves give h(x,t) for t = 0:2:8 with vertical offsets of',...
           ' t/4.'])
end
hold off

%=========================================================================
% Subfunctions

    function alpha = getalpha(t,epsilon,beta)
    % Compute alpha by iteration.
        ao = 2*t + beta;
        for i = 1:10
            alpha = 2*t + beta ...
                    - 2*epsilon*( -0.375*(ao-beta).*sin(beta)...
                    + 0.125*(cos(beta) - cos(ao)) )...
                 - (2/128)*epsilon^2*( (ao-beta).*(14 - 15*cos(2*beta))...
                    + (-13*sin(beta).*cos(beta) + sin(ao).*cos(ao) ...
                    + 12*sin(beta).*cos(ao)) );

            if norm(alpha - ao,inf) <= 0.001
                return
            else
                ao = alpha;
            end
        end
        error(['Unable to compute alpha for t = ',num2str(t),'.'])
    
    % end function getalpha

% end function strainedA