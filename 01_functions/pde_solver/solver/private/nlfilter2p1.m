function v = nlfilter2p1(u)
    % Nonlinear filter Algorithm 2.1 of Engquist et alia.

    v = u;
    for j = 2:length(v)-1
        Dp = v(j+1) - v(j);
        Dm = v(j) - v(j-1);
        if Dp*Dm < 0 
            if abs(Dp) > abs(Dm)
                corr = sign(Dp)*min(abs(Dm),abs(Dp)/2);
                v(j) = v(j) + corr;
                v(j+1) = v(j+1) - corr;
            else
                corr = sign(Dp)*min(abs(Dp),abs(Dm)/2);
                v(j) = v(j) + corr;
                v(j-1) = v(j-1) - corr;            
            end
        end
    end
        
% end function nlfilter2p1