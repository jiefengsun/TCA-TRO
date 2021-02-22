function [l_t, l_star, r_star, alpha_star,r_t, N , alpha_min] = TCA_geo( mw)
% 6/05/2020 creat this script
% Cannot pass varialbe to dynamic work space (cannot load it using run), 
% Wrap it as a function and pass varialbes as arguments. - It seems
% clearer. 
% 12/7/2020 modified from TCA_parameters to only provide geometry parameters
% 12/18/2020 add the stress output
% 2/12/2020 revised output, removed theta_bar_star
  
% Independent parameters
    l_t =  175e-3; %
    r_t = 0.21e-3;
    r_m = 0.406e-3;
    alpha_m = 22.42/180*pi;
    l_m = 66.75e-3; 
    N = sqrt(l_t^2-l_m^2)/(2*pi*r_m); % The measured value is 63; 
    
    
% creeped length table
    Mw = [0, 2, 30, 60]; 
    Alpha_star_deg = [15, 16, 20, 22];
    alpha_star_deg = interp1(Mw, Alpha_star_deg, mw); 


    alpha_star = alpha_star_deg/180*pi; % rad corresponding to the instant equlibrium when weigth is removed. 
    r_star = r_m*cos(alpha_star) /cos(alpha_m); 
     
%  Fixed parameter
   
    l_star = l_t*sin(alpha_star); %    
     alpha_min = 12/180*pi; % the pitch angle that some coils starts to contact   

    
end