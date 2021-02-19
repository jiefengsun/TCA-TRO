function out = ad(xi)
    % this is Adjoint representation of se(3)
    % xi is the input 6x1 vector, the out put will be homogeneous transformation matrix

    out = [skew(xi(1:3)),zeros(3);skew(xi(4:6)),skew(xi(1:3))];
end