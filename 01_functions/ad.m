function adV= ad(V)
    % this is Adjoint representation of se(3)
    % take: V that is 6x1 vector, 
    % output:  the out put will be a 6x6 matrix
    
    % Used to calculate the Lie bracket [V1, V2] = [adV1]V2
    % Example Input:
    %{ 
      clear;clc;
      V = [1; 2; 3; 4; 5; 6];
      adV = ad(V)
    %}
    % Output:
    % adV =
    %     0    -3     2     0     0     0
    %     3     0    -1     0     0     0
    %    -2     1     0     0     0     0
    %     0    -6     5     0    -3     2
    %     6     0    -4     3     0    -1
    %    -5     4     0    -2     1     0
   
    skew_m = skew(V(1:3));
    adV = [skew_m ,zeros(3);skew(V(4:6)),skew_m ];
end