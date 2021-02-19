function out = skew(u)
    % to calculate the skew syemmetric matrix, input u is a vector

    out =  [0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
end