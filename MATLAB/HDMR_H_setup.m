function [N,d,graphics,maxorder,maxiter,bf1,m,K,R,method,alfa, ...
    vartol,lambda,refit] = HDMR_H_setup(X,y,user_options)
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
% Check/setup of the input arguments of HDMR_H function                   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [N,d,graphics,maxorder,maxiter,bf1,m,K,R,method,alfa,vartol, ...      %
%       lambda,refit] = HDMR_H_setup(X,y,user_options)                    %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Determine the size of matrix X
[N,d] = size(X); 
if ( d == 1 )
    error(['HDMR_H ERROR: Matrix X contains only a single column: ' ...
        'No point to do sensitivity analysis when d = 1\n']); 
end
if ( N < 300 )
    error(['HDMR_H ERROR: Number of samples, N, of N x d matrix X ' ...
        'is insufficient\n']); 
end
[a,b] = size(y);
if ( a ~= N )
   error(['HDMR_H ERROR: Dimension mismatch. The number of rows, N ' ...
       'of Y should match number of rows of X\n']); 
end 
if ( b ~= 1 )
   error(['HDMR_H ERROR: Y should be a N x 1 vector with one ' ...
       'simulated output for each N parameter vectors\n']); 
end 

% Define default options
def_options = struct('graphics',1,'maxorder',3,'maxiter',100, ...
    'bf1',1,'m',2,'K',100,'R',floor(N/2),'method',1,'alfa',0.01, ...
    'lambda',0.1,'vartol',1e-5,'refit',1);

% Open file for settings.txt
fid = fopen('HDMR_H_settings.txt','w');
fprintf(fid,'|-------------------------------------------------------------------------|\n');
fprintf(fid,'|                                                                         |\n');
fprintf(fid,'|   HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR        HHH     HHH  |\n');
fprintf(fid,'|   HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR       HHH     HHH  |\n');
fprintf(fid,'|   HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR      HHH     HHH  |\n');
fprintf(fid,'|   HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR      HHH     HHH  |\n');
fprintf(fid,'|   HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR   --- HHHHHHHHHHH  |\n');
fprintf(fid,'|   HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR    --- HHHHHHHHHHH  |\n');
fprintf(fid,'|   HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR       HHH     HHH  |\n');
fprintf(fid,'|   HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR       HHH     HHH  |\n');
fprintf(fid,'|   HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR      HHH     HHH  |\n');
fprintf(fid,'|   HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR     HHH     HHH  |\n');
fprintf(fid,'|                                                                         |\n');
fprintf(fid,'|  SSSSSSS EEEEEEEE TTTTTTTT TTTTTTTT IIIIIIII NN    NN  GGGGGG   SSSSSSS |\n');
fprintf(fid,'| SSSSSSS  EEEEEEEE TTTTTTTT TTTTTTTT  IIIIII  NNN   NN GGGGGG   SSSSSSS  |\n');
fprintf(fid,'| SS       EE       TT TT TT TT TT TT    II    NNNN  NN GG       SS       |\n');
fprintf(fid,'| SSSSSS   EEEEE    T  TT  T T  TT  T    II    NN NN NN GG  GGGG SSSSSS   |\n');
fprintf(fid,'| SSSSSSS  EEEEE       TT       TT       II    NN NN NN GG   GGG SSSSSSS  |\n');
fprintf(fid,'|       SS EE          TT       TT       II    NN  NNNN GG    GG       SS |\n');
fprintf(fid,'|  SSSSSSS EEEEEEEE    TT       TT     IIIIII  NN   NNN GGGGGGGG  SSSSSSS |\n');
fprintf(fid,'| SSSSSSS  EEEEEEEE    TT       TT    IIIIIIII NN    NN  GGGGGG  SSSSSSS  |\n');
fprintf(fid,'|                                                                         |\n');
fprintf(fid,'|-------------------------------------------------------------------------|\n');

% Now unpack structure user_options
if isfield(user_options,'graphics')
    graphics = user_options.graphics;
    if ~any(ismember(graphics,[0 1]))
        error(['HDMR_H ERROR: Field "graphics" of options should take ' ...
            'on the value of 0 or 1\n']);
    end
else
    fprintf(['HDMR_H DEFAULT: Field "graphics" of options not ' ...
        'specified. Default: graphics = 1 assumed\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "graphics" of options not ' ...
        'specified. Default: graphics = 1 assumed\n']);
    graphics = def_options.graphics;
end

if isfield(user_options,'maxorder')
    maxorder = user_options.maxorder;
    if ~any(ismember(maxorder,[1 2 3]))
        error(['HDMR_H ERROR: Field "maxorder" of options should be ' ...
            'an integer with values of 1, 2 or 3\n']);
    end
else
    fprintf(['HDMR_H DEFAULT: Field "maxorder" of options not ' ...
        'specified. Default: maxorder = 3 used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "maxorder" of options not ' ...
        'specified. Default: maxorder = 3 used\n']);
    maxorder = def_options.maxorder;
end
% Important next check for maxorder - as maxorder relates to d
if ( d == 2 )
    if maxorder > 2
        warning(['HDMR_H WARNING: Field "maxorder" of options set ' ...
            'to 2 as d = 2 (X has two columns) \n']);
        maxorder = 2;
    end
end
if isfield(user_options,'maxiter')
    maxiter = user_options.maxiter;
    if ~any(ismember(maxiter,1:1000))
        error(['HDMR_H ERROR: Field "maxiter" of options should be ' ...
            'an integer between 1 to 1000\n']);
    end
else
    fprintf(['HDMR_H DEFAULT: Field "maxiter" of options not ' ...
        'specified. Default: maxiter = 100 used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "maxiter" of options not ' ...
        'specified. Default: maxiter = 100 used\n']);
    maxiter = def_options.maxiter;
end
if isfield(user_options,'bf1')
    bf1 = user_options.bf1;
    if ~ismember(bf1,[0 1])
        error(['HDMR_H ERROR: Field "bf1" of options should be an ' ...
            'integer with value 0 or 1\n']);
    end
else
    fprintf(['HDMR_H DEFAULT: Field "bf1" of options not ' ...
        'specified. Default: bf1 = 1 is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "bf1" of options not ' ...
        'specified. Default: bf1 = 1 is used\n']);
    bf1 = def_options.bf1;
end
if isfield(user_options,'m')
    m = user_options.m;
    if ~any(ismember(m,1:10))
        error(['HDMR_H ERROR: Field "m" of options should be an ' ...
            'integer between 1 to 10\n']);
    end    
else
    fprintf(['HDMR_H DEFAULT: Field "m" of options not ' ...
        'specified. Default: m = 2 is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "m" of options not ' ...
        'specified. Default: m = 2 is used\n']);
    m = def_options.m;
end
if isfield(user_options,'K')
    K = user_options.K;
    if ~any(ismember(K,1:500))
        error(['HDMR_H ERROR: Field "K" of options should be an ' ...
            'integer between 1 to 500\n']);
    end        
else
    fprintf(['HDMR_H DEFAULT: Field "K" of options not ' ...
        'specified. Default: K = 100 is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "K" of options not ' ...
        'specified. Default: K = 100 is used\n']);
    K = def_options.K;
end
if isfield(user_options,'R')
    R = user_options.R;
    if ~any(ismember(R,100:N))
        error(['HDMR_H ERROR: Field "R" of options should be an ' ...
            'integer between 100 and N, the number of rows matrix X\n']);
    end            
else
    fprintf(['HDMR_H DEFAULT: Field "R" of options not ' ...
        'specified. Default: R = floor(N/2) is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "R" of options not ' ...
        'specified. Default: R = floor(N/2) is used\n']);
    R = def_options.R;
end
if isfield(user_options,'method')
    method = user_options.method;
    if ~any(ismember(method,[1 2]))
        error(['HDMR_H ERROR: Field "method" of options should take ' ...
            'on the value of 1 (forward selection) or ' ...
            '2 (backward elimination)\n']);
    end    
else
    fprintf(['HDMR_H DEFAULT: Field "method" of options not ' ...
        'specified. Default: method = 1 (forward selection) is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "method" of options not ' ...
        'specified. Default: method = 1 (forward selection) is used\n']);
    method = def_options.method;
end
if isfield(user_options,'alfa')
    alfa = user_options.alfa;
    if ischar(alfa)
        error(['HDMR_H ERROR: Field "alfa" (confidence interval) of ' ...
            'options should not be a string but a numerical value\n']);
    end 
    if alfa > 0.5
        error(['HDMR_H ERROR: Field "alfa" (confidence interval) of ' ...
            'options should not be larger than 0.5. Default: 0.01\n']);
    end    
    if alfa < eps
        error(['HDMR_H ERROR: Field "alfa" (confidence interval) of ' ...
            'options should at least be equal to eps. Default: 0.01\n']);
    end    
else
    fprintf(['HDMR_H DEFAULT: Field "alfa" of options not ' ...
        'specified. Default: alfa = 0.01 is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "alfa" of options not ' ...
        'specified. Default: alfa = 0.01 is used\n']);
    alfa = def_options.alfa;
end
if isfield(user_options,'lambda')
    lambda = user_options.lambda;
    if ischar(lambda)
        error(['HDMR_H ERROR: Field "lambda" (regularization term) of ' ...
            'options should not be a string but a numerical value\n']);
    end 
    if lambda > 10
        fprintf(['HDMR_H WARNING: Field "lambda" of options set ' ...
            'rather large. Default: lambda = 0.1\n']);
        fprintf(fid,['HDMR_H WARNING: Field "lambda" of options ' ...
            'set rather large. Default: lambda = 0.1\n']);
    end    
    if lambda < 0
        error(['HDMR_H ERROR: Field "lambda" (regularization term) ' ...
            'of options cannot be smaller than zero. Default: 0.1\n']);
    end
else
    fprintf(['HDMR_H DEFAULT: Field "lambda" of options not ' ...
        'specified. Default: lambda = 0.1 is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "lambda" of options not ' ...
        'specified. Default: lambda = 0.1 is used\n']);
    lambda = def_options.lambda;
end
if isfield(user_options,'vartol')
    vartol = user_options.vartol;
    if ischar(vartol)
        error(['HDMR_H ERROR: Field "vartol" (convergence threshold) ' ...
            'of options should not be a string but a numerical value\n']);
    end
    if vartol < 0
        error(['HDMR_H ERROR: Field "vartol" (convergence threshold) ' ...
            'of options cannot be negative. Default: 1e-5\n']);
    end    
    if vartol < 1e-8
        error(['HDMR_H ERROR: Field "vartol" (convergence threshold) ' ...
            'of options set rather small. Default: 1e-5\n']);
    end    
    if vartol > 0.1
        error(['HDMR_H ERROR: Field "vartol" (convergence threshold) ' ...
            'of options set rather large. Default: 1e-5\n']);
    end    
else
    fprintf(['HDMR_H DEFAULT: Field "vartol" of options not ' ...
        'specified. Default: vartol = 1e-3 is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "vartol" of options not ' ...
        'specified. Default: vartol = 1e-3 is used\n']);
    vartol = def_options.vartol;
end
if isfield(user_options,'refit')
    refit = user_options.refit;
    if ~any(ismember(refit,[0 1]))
        error(['HDMR_H ERROR: Field "refit" of options should take on ' ...
            'the value of 0 or 1\n']);
    end
else
    fprintf(['HDMR_H DEFAULT: Field "refit" of options not ' ...
        'specified. Default: refit = 1 is used\n']);
    fprintf(fid,['HDMR_H DEFAULT: Field "refit" of options not ' ...
        'specified. Default: refit = 1 is used\n']);
    refit = def_options.refit;
end

fprintf(fid,'\n');
fprintf(fid,'          |===================================|\n');
fprintf(fid,'          |  field of options      value      |\n');
fprintf(fid,'          |-----------------------------------|\n');
fprintf(fid,'          |    %-12s \t %8d     |\n','graphics   ',graphics);
fprintf(fid,'          |    %-12s \t %8d     |\n','maxorder   ',maxorder);
fprintf(fid,'          |    %-12s \t %8d     |\n','maxiter    ',maxiter);
fprintf(fid,'          |    %-12s \t %8d     |\n','bf1        ',bf1);
fprintf(fid,'          |    %-12s \t %8d     |\n','method     ',method);
fprintf(fid,'          |    %-12s \t %8d     |\n','m          ',m);
fprintf(fid,'          |    %-12s \t %8d     |\n','K          ',K);
fprintf(fid,'          |    %-12s \t %8d     |\n','R          ',R);
fprintf(fid,'          |    %-12s \t %8f     |\n','alfa       ',alfa);
fprintf(fid,'          |    %-12s \t %8.2e     |\n','lambda     ',lambda);
fprintf(fid,'          |    %-12s \t %8.2e     |\n','vartol     ',vartol);
fprintf(fid,'          |===================================|\n');
fprintf(fid,'\n');
% Now close file
fclose(fid);

end
