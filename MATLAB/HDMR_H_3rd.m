function [y_ijk,C3,it] = HDMR_H_3rd(B1,y_res,R,gamma,n3,m3,c3,lambda)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR        HHH     HHH   %
%  HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR       HHH     HHH   %
%  HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR      HHH     HHH   %
%  HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR      HHH     HHH   %
%  HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR   --- HHHHHHHHHHH   %
%  HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR    --- HHHHHHHHHHH   %
%  HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR       HHH     HHH   %
%  HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR       HHH     HHH   %
%  HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR      HHH     HHH   %
%  HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR     HHH     HHH   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Third order - individual estimation and backfitting                     %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [y_ijk,C3,it] = HDMR_H_3rd(B1,y_res,R,gamma,n3,m3,c3,lambda)          %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

B3 = zeros(R,m3);       % Initialize B2 matrix
C3 = zeros(m3,n3);      % Initialize coefficients
T3 = zeros(m3,R);       % Initialize T(emporary) matrix - 3rd
y_ijk = zeros(R,n3);    % Initialize 3rd order contributions
it = 0;                 % Initialize iteration counter

% Third order individual estimation
for j = 1:n3
    % First compute the B-spline coefficients
    for z = 1:m3
        B3(1:R,z) = B1(:,gamma(z,1),c3(j,1)) .* ...
            B1(:,gamma(z,2),c3(j,2)) .* B1(:,gamma(z,3),c3(j,3));
    end
    % Regularized least squares inversion ( minimize || C3 ||_2 )
    % C3(:,k) = ( B3(1:R,:,k)'*B3(1:R,:,k) + ...
    %               lambda * eye(m3) ) \ B3(1:R,:,k)' * y_res;
    B33 = B3' * B3;
    if det(B33) == 0
        % Ill-defined --> C3(:,k) = 0;
    else
        T3 = ( B33 + lambda * eye(m3) ) \ B3';
    end
    % Derive coefficients
    C3(1:m3,j) = T3 * y_res;
    % Derive jth term of 3rd order contribution of emulator
    y_ijk(1:R,j) = B3 * C3(1:m3,j);
end

end
