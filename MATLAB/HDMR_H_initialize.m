function [Xy,Em,SA,RT,Y_em,beta,gamma,m1,m2,m3,j1,j2,j3,it2, ...
    it3,itr] = HDMR_H_initialize(X,y,N,d,K,R,m,maxorder)
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
% Initialize main variables used by HDMR                                  %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [Xy,Em,SA,RT,Y_em,beta,gamma,m1,m2,m3,j1,j2,j3,it2,it3,itr] = ...     %
%       HDMR_H_initialize(X,y,N,d,K,R,m,maxorder)                         %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Random seed (legacy: randn('state',sum(100*clock)); )
rng(1+round(100*rand),'twister');

% STRUCTURE XY: Define content
if K == 1
    Xy = struct('X_n',nan(N,d),'minX',min(X,[],1),'maxX',max(X,[],1), ...
        'Y',y,'R',R,'id',(1:R)');
else
    % Now setup the boostrap matrix with selection matrix, id, for samples
    [~,id] = sort(rand(N,K));
    % Store in structure XY
    Xy = struct('X_n',nan(N,d),'minX',min(X,[],1),'maxX',max(X,[],1), ...
        'Y',y,'R',R,'id',id(1:R,1:K));
end
% Compute normalized X-values
Xy.X_n = (X(1:N,1:d) - repmat(Xy.minX(1:d),N,1)) ./ ...
    repmat(Xy.maxX(1:d) - Xy.minX(1:d),N,1);

% STRUCTURE Em: Important variables
[n2,n3] = deal(0); [c2,c3] = deal([]);
% determine all combinations of parameters (coefficients) for each order
c1 = 1:d; n1 = d; % = size(c1,1);
% now return based on maxorder used
if maxorder > 1
    c2 = nchoosek(1:d,2); n2 = size(c2,1);
end
if maxorder == 3
    c3 = nchoosek(1:d,3); n3 = size(c3,1);
end
% calulate total number of coefficients
n = n1 + n2 + n3;

% Initialize m1, m2 and m3 - # coefficients first, second, third order
m1 = m + 3; m2 = m1^2; m3 = m1^3;

% STRUCTURE Em: Initialization
switch maxorder
    case 1
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),'c1',c1,'n1',d, ...
            'c2',c2,'n2',n2,'c3',c3,'n3',n3,'n',n,'maxorder',maxorder, ...
            'select',nan(n,K),'C1',nan(m1,n1,K),'C2',nan(1,1,K), ...
            'C3',nan(1,1,K),'B1',zeros(N,m1,n1),'iter',nan(4,K));
    case 2
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),'c1',c1,'n1',n1, ...
            'c2',c2,'n2',n2,'c3',c3,'n3',n3,'n',n,'maxorder',maxorder, ...
            'select',nan(n,K),'C1',nan(m1,n1,K),'C2',nan(m2,n2,K), ...
            'C3',nan(1,1,K),'B1',zeros(N,m1,n1),'iter',nan(4,K));
    case 3
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),'c1',c1,'n1',n1, ...
            'c2',c2,'n2',n2,'c3',c3,'n3',n3,'n',n,'maxorder',maxorder, ...
            'select',nan(n,K),'C1',nan(m1,n1,K),'C2',nan(m2,n2,K), ...
            'C3',nan(m3,n3,K),'B1',zeros(N,m1,n1),'iter',nan(4,K));
end

% Now compute B-spline values for all N samples of X_n
Em.B1 = B_spline(Xy.X_n,N,d,m); 

% STRUCTURE SA: Sensitivity analysis and analysis of variance decomposition
SA = struct('S',nan(Em.n,K),'Sa',nan(Em.n,K),'Sb',nan(Em.n,K), ...
    'ST',nan(d,1),'V_em',nan(Em.n,K),'V_y',nan(1,K));
% Return runtime
RT = zeros(1,K);
% Initialize various terms
Y_em = nan(R,Em.n);
% Initialize beta and gamma as 2nd and 3rd order combinations
beta = permn(1:m1,2); gamma = permn(1:m1,3);
% Initialize index of first, second and third order terms (columns of Y_bf)
j1 = (1:n1)'; j2 = (n1+1:n1+n2)'; j3 = (n1+n2+1:n)';
% Now initialize number of iterations first, second and third order
[it2,it3,itr] = deal(0);

fprintf('\n');
fprintf('  ------------------------------------------------------------------------------  \n');
fprintf('       HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR        HHH     HHH       \n');
fprintf('       HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR       HHH     HHH       \n');
fprintf('       HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR      HHH     HHH       \n');
fprintf('       HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR      HHH     HHH       \n');
fprintf('       HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR   --- HHHHHHHHHHH       \n');
fprintf('       HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR    --- HHHHHHHHHHH       \n');
fprintf('       HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR       HHH     HHH       \n');
fprintf('       HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR       HHH     HHH       \n');
fprintf('       HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR      HHH     HHH       \n');
fprintf('       HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR     HHH     HHH       \n');
fprintf('  ------------------------------------------------------------------------------  \n');
fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
fprintf('\n'); 

end
