function [t, T] = thermal_model(r_t, l_t, dt, input_data )

% the thermal model that will be used for the dynamics of TCAs. 
% The input parameters: [r_t, l_t, dt, input_data] are respectively the 
% radius of a twisted fiber, the length of a twisted fiber, the time step, 
% and input power 

% initialize parameters
rho = 1127; % kg/m^3, this is corrected. 
cp = 2000; % J/(K kg), silver 2400, nylon 1700, so the silver nylon will be at the middle. 
At = r_t^2*pi; 
V = l_t*At*rho; % m^3
mt = rho*V/1000; %kg
h = 23; %% used 23 for dynamics, used 32 for varying load
As = 3*2*pi*r_t*l_t; % surface is enlarged. 
T0  = 25; 
% sigma = 5.67e-8; %Boltzmann constantW/m2 K4
% epsilon = 0.5; % emissivity.
Pin_data = readmatrix(input_data); %'time_Pin.txt'
%  interpolate the input power and turn the unit to seconds. 
tspan = 0:dt:max(Pin_data(:,1))/1000; 
% plot(0:0.1:4,  interp1(Pin_data(:, 1)/1000, Pin_data(:, 2), tspan))
[t, T]= ode45(@myode, tspan, T0); 
plot(t, T)
drawnow

    function dTdt = myode(t, T)          
         dTdt=  1/(mt*cp)*(-h*As*(T-T0)+ interp1(Pin_data(:, 1)/1000, Pin_data(:, 2), t));
    end


end





