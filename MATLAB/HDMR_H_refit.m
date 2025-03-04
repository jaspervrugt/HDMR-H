function [Y_em,C1,C2,C3,iter] = HDMR_H_refit(B1,T1,C1,C2,C3,y_res, ...
    Y_em,R,beta,gamma,n1,n2,n,c2,c3,m1,m2,m3,j1,j2,j3,ind,lambda, ...
    vartol,maxiter,maxorder)
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
% Set to zero insignificant components in coefficient matrix + Y_em and   %
% redo backfitting                                                        %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [Y_em,C1,C2,C3,iter] = HDMR_H_refit(B1,T1,C1,C2,C3,y_res,Y_em, ...    %
%       R,beta,gamma,n1,n2,n,c2,c3,m1,m2,m3,j1,j2,j3,ind,lambda, ...      %
%       vartol,maxiter,maxorder)                                          %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Initialize important variables used for convergence analysis backfitting
varb_old = nan(1,n); varb_max = 1; iter = 0; i2 = []; i3 = [];

i1 = find(ind(j1) == 0);                    % indx insignfcnt 1st ord terms
C1(1:m1,i1) = 0; Y_em(1:R,j1(i1)) = 0;      % set id columns to zero
varb_old(1,j1) = sum(C1.^2,1);              % compute old sum coefficients

if maxorder > 1
    i2 = find(ind(j2) == 0);                % indx insignfcnt 2nd ord terms
    C2(1:m2,i2) = 0; Y_em(1:R,j2(i2)) = 0;  % set id columns to zero
    varb_old(1,j2) = sum(C2.^2,1);          % compute old sum coefficients
    % ii2 = n2 - numel(i2);                   % # significant 2nd order terms
    % B2 = nan(R,m2,ii2); T2 = nan(m2,R,ii2); % initialize B2 and T2 matrices
end

if maxorder > 2
    i3 = find(ind(j3) == 0);                % indx insignfcnt 3rd ord terms
    C3(1:m3,i3) = 0; Y_em(1:R,j3(i3)) = 0;  % set id columns to zero
    varb_old(1,j3) = sum(C3.^2,1);          % compute old sum coefficients
    % ii3 = n3 - numel(i3);                   % # significant 3rd order terms
    % B3 = nan(R,m3,ii3); T3 = nan(m3,R,ii3); % initialize B3 and T3 matrices    
end

% Now collect all coefficients with insignificant sensitivity
ii = [ i1 ; j2(i2) ; j3(i3) ];
% Now determine index of all significant terms
idx = 1:n; idx(ii) = [];
% Initialize B2 and B3 matrices
B2 = zeros(R,m2); B3 = zeros(R,m3);

% Backfitting of all orders combined
while ( varb_max > vartol ) && ( iter < maxiter )
    varb_new = zeros(1,n); ct_2 = 1; ct_3 = 1;
    for i = 1:numel(idx)
        % Get column and index of all other non-zero columns
        j = idx(i); ii = idx; ii(i) = [];
        % Now remove all non-zero columns from Y_em except column j
        Y_new(1:R,1) = y_res(1:R,1) - sum(Y_em(1:R,ii),2);
        if ismember(j,j1)
            C1(1:m1,j) = T1(1:m1,1:R,j) * Y_new(1:R,1);
            Y_em(1:R,j) = B1(1:R,1:m1,j) * C1(1:m1,j);
            varb_new(j) = sum(C1(1:m1,j).^2);
        elseif ismember(j,j2)
            z = j - n1;
            % Store T2 of significant 2nd order terms at first iteration
            if (iter == 0)
                for q = 1:m2
                    B2(1:R,q) = B1(1:R,beta(q,1),c2(z,1)) .* ...
                        B1(1:R,beta(q,2),c2(z,2));
                end
                T2(1:m2,1:R) = ( B2(1:R,1:m2)' * ...
                    B2(1:R,1:m2) + lambda * eye(m2) ) ...
                    \ B2(1:R,1:m2)';
            end
            C2(1:m2,z) = T2(1:m2,1:R) * Y_new(1:R,1);
            Y_em(1:R,j) = B2(1:R,1:m2) * C2(1:m2,z);
            varb_new(j) = sum(C2(1:m2,z).^2);
            ct_2 = ct_2 + 1;
        else
            z = j - n1 - n2;
            % Store T3 of significant 3rd order terms at first iteration
            if (iter == 0)
                for q = 1:m3
                    B3(1:R,q) = B1(1:R,gamma(q,1),c3(z,1)) .* ...
                        B1(1:R,gamma(q,2),c3(z,2)) .* ...
                        B1(1:R,gamma(q,3),c3(z,3));
                end
                T3(1:m3,1:R) = ( B3(1:R,1:m3)' * ...
                    B3(1:R,1:m3) + lambda * eye(m3) ) ...
                    \ B3(1:R,1:m3)';
            end
            C3(1:m3,z) = T3(1:m3,1:R) * Y_new(1:R,1);
            Y_em(1:R,j) = B3(1:R,1:m3) * C3(1:m3,z);
            varb_new(j) = sum(C3(1:m3,z).^2);
            ct_3 = ct_3 + 1;
        end
    end
    varb_max = max(abs(varb_new - varb_old));
    varb_old = varb_new;
    iter = iter + 1;
end

end
