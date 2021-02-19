function out = unskew(m)
% undo the skew
% convert a skew symmetric matrix to a vector
    out = [m(3,2);m(1,3);m(2,1)];
end