function out = unse(m)
    out = [unskew(m(1:3,1:3));m(1:3,4)];
end