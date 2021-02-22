function out = skew(u)
    % Takes a 3x1-vector  
    % Returns the skew symmetric matrix in so(3).
    % Example Input:
    %{
      clear;clc;
      u = [1; 2; 3];
      skew = VecToso3(u)
    %}
    % Output:
    % skew =
    %     0    -3     2
    %     3     0    -1
    %    -2     1     0

    out =  [0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
end