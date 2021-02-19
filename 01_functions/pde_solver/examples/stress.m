function stress
% Problem of section 9, "Waves on a Traveling Threadline",
% somewhat modified. W = (m,V)^T

t = 0;
Npoints = 1000;
x = linspace(-10,10,Npoints);
W = ones(2,Npoints);
W(2,:) = exp(-abs(x));

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
if meth == 4
    form = 3;
    pdefun = @cl;
else
    % FORM = 2 has not been coded for these PDEs.
    mform = menu('Specify the form of the PDEs','1','3');
    if mform == 1
        form = 1;
        pdefun = @pdes;
    else
        form = 3;
        pdefun = @cl;
    end
end

periodic = true;
sol = setup(form,pdefun,t,x,W,method,periodic);

close all
m = sol.u(1,:);
V = sol.u(2,:);
plot(x,m,'r',x,V,'b')
legend('m(x)','V(x)',2)
xlabel(['Solution at t = ',num2str(t),'.']);
axis([-10 10 -0.1 2])
pause

howfar = 1;
for m = 1:5    
    sol = hpde(sol,howfar,@timestep);
    t = sol.t;
    m = sol.u(1,:);
    V = sol.u(2,:);
    plot(x,m,'r',x,V,'b')
    legend('m(x)','V(x)',2)
    title(['Computed with ',method,'.'])
    xlabel(['Solution at t = ',num2str(t),'.']);
    axis([-10 10 -0.1 2])
    pause
end

%=========================================================================
% Subfunctions

    function F = pdes(t,x,W,W_x)
        % W = (m,V)^T
        F = zeros(size(W));
        F(1,:) = - [W(1,:).*W_x(2,:) + W(2,:).*W_x(1,:)];
        F(2,:) = - [W_x(1,:)./W(1,:).^3 + W(2,:).*W_x(2,:)];
    % end function pdes
    
    function F = cl(W)
        % W = (m,V)^T        
        F = [-W(1,:).*W(2,:); (1./W(1,:).^2 - W(2,:).^2)/2];
    % end function cl

    function dt = timestep(dx,t,x,W)
        % W = (m,V)^T
        mumax1 = max(abs(W(2,:) + 1./W(1,:)));
        mumax2 = max(abs(W(2,:) - 1./W(1,:)));
        dt = 0.9*dx/max(mumax1,mumax2);
    % end function timestep

% end function stress