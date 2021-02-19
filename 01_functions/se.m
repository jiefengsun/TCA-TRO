function out = se(x)
    %This is the the hat operator. 
    out = [skew(x(1:3)),x(4:6);0,0,0,0];
end