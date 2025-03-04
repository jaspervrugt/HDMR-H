function [y_res,y_ij,C2,it] = HDMR_H_2nd(B1,y_res,R,beta,n2,m2,c2,lambda)
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
% Second order - individual estimation and backfitting                    %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [y_res,y_ij,C2,it] = HDMR_H_2nd(B1,y_res,R,beta,n2,m2,c2,lambda)      %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

B2 = zeros(R,m2);       % Initialize B2 matrix
C2 = zeros(m2,n2);      % Initialize coefficients
T2 = zeros(m2,R);       % Initialize T(emporary) matrix - 2nd
y_ij = zeros(R,n2);     % Initialize 2nd order contributions
it = 0;                 % Initialize iteration counter

% Second order individual estimation
for j = 1:n2
    % First compute the B-spline coefficients
    for z = 1:m2
        B2(1:R,z) = B1(1:R,beta(z,1),c2(j,1)) .* B1(1:R,beta(z,2),c2(j,2));
    end
    % Regularized least squares inversion ( minimize || C2 ||_2 )
    % C2(:,k) = ( B2(1:R,:,k)'*B2(1:R,:,k) + ...
    %               lambda * eye(m2) ) \ B2(1:R,:,k)' * y_res;
    B22 = B2'*B2;
    if det(B22) == 0
        % Ill-defined --> C2(:,k) = 0;
    else
        T2 = ( B22 + lambda * eye(m2) ) \ B2';
    end
    % Derive coefficients
    C2(1:m2,j) = T2 * y_res;
    % Derive jth term of 2nd order contribution of emulator    
    y_ij(1:R,j) = B2 * C2(1:m2,j);
end

end
