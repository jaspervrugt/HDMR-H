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
%                                                                         %
%  High-Dimensional Model Representation (HDMR) for d > 20 using B-spline %
%  functions for variance-based global sensitivity analysis (GSA) with    %
%  correlated and uncorrelated inputs. This function uses as input a      %
%  N x d matrix of N different d-vectors of model inputs (parameters) and %
%  a N x 1 vector of corresponding model outputs and returns to the user  %
%  each factor's first, second, and third order sensitivity coefficient   %
%  (separated in total, structural and correlative contributions), an     %
%  estimate of their 95% confidence intervals (from bootstrap method)     %
%  and the coefficients of the significant B-spline basis functions that  %
%  govern output, Y (determined by an F-test of the error residuals of    %
%  the HDMR model (emulator) with/without a given first, second and/or    %
%  third order B-spline). These coefficients define an emulator that can  %
%  be used to predict the output, Y, of the original (CPU-intensive?)     %
%  model for any d-vector of model inputs. For uncorrelated model inputs  %
%  (columns of X are independent), the HDMR sensitivity indices reduce    %
%  to a single index (= structural contribution), consistent with their   %
%  values derived from commonly used variance-based GSA methods.          %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [S,Ss,Fx,Em,Xy,RT] = HDMR_H(X,y);                                     %
%   [S,Ss,Fx,Em,Xy,RT] = HDMR_H(X,y,options);                             %
%                                                                         %
%  INPUT ARGUMENTS                                                        %
%   X           Nxd matrix: N vectors of d parameters                     %
%   y           Nx1 vector: single model output for each row of X         %
%   options     (optional) structure: Fields specify HDMR variables       %
%    .graphics  integer [0,1]: graphical output?             (def: 1)     %
%    .maxorder  integer [1-3]: max order emulator            (def: 3)     %
%    .maxiter   integer [1-1000]: max iterations backfitting (def: 100)   %
%    .bf1       integer [0,1]: 1st order, ind./backfitting   (def: 1)     %
%    .method    integer [1,2]: 1=forw. sel.; 2=backw. elim.  (def: 1)     %
%    .m         integer [2-10]: # B-spline intervals         (def: 2)     %
%    .K         integer [1-500] # bootstrap iterations       (def: 100)   %
%    .R         integer [100-N/2] # bootstrap samples        (def: N/2)   %
%    .alfa      real [0.-0.1]: significance level F-test     (def: 0.01)  %
%    .lambda    real [0-inf]: regularization coefficient     (def: 0.10)  %
%    .vartol    real (0-1]: tolerance backfitting            (def: 1e-5)  %
%                                                                         %
%   def_options = struct('graphics',1,'maxorder',3,'maxiter','100',...    %
%                   'bf1',1,'method',1,'m',2,'K',100,'R',N/2,...          %
%                   'alfa',0.01,'lambda',0.1,'vartol',1e-5);              %
%                                                                         %
%  OUTPUT ARGUMENTS                                                       %
%   S          Cell array: structural, correlative and total sensitivity  %
%              each component function -> 1st, 2nd and 3rd order effects  %
%              ( = Table 2 on screen and printed in "HDMR_results.txt" )  %
%   Ss         Cell array: As "S" but lists only significant component    %
%              functions determined via model selection using a F-test    %
%              ( = Table 3 on screen and printed in "HDMR_results.txt" )  %
%   Fx         Cell array: Tabulates emulator properties and performance  %
%              on training data set for each of the K bootstrap trials    %
%              ( = Table 1 on screen and printed in "HDMR_results.txt" )  %
%   Em         Structure array: Fields with input and output K emulators  %
%    .B1       Nx(m+3)xn1 matrix: B-spline function evaluated at X        %
%    .C1       dxn1xK matrix: 1st order coefficients (via backfitting)    %
%    .C2       d2xn2xK matrix: 2nd order coefficients (via backfiting)    %
%    .C3       d^3xn3xK matrix: 3rd order coefficients (via backfitting)  %
%    .m        scalar: number of spline intervals ( = m + 3)              %
%    .Y_e      RxK matrix: Emulator predictions of K training data sets   %
%    .RMSE     1xK vector: RMSE of emulator residuals of K training sets  %
%    .c1       n1x1 vector: indices of 1st order terms ( = 1:d )          %
%    .c2       n2x2 matrix: indices of 2nd order combinations             %
%    .c3       n3x3 matrix: indices of 3rd order combinations             %
%    .f0       1xK vector: mean Y of each of the K bootstrap trials       %
%    .n        scalar: ( = n1+n2+n3 ) max. number component functions     %
%    .n1       scalar: ( = d) with total number of 1st order terms        %
%    .n2       scalar: total number of 2nd order componen functions       %
%    .n3       scalar: total number of 3rd order component functions      %
%    .maxorder scalar: maximum order of HDMR emulator                     %
%    .nterms   1xK vector: # significant terms of each of K emulators     %
%    .p0       1xK vector: # parameters of each of the K emulators        %
%    .select   nxK matrix: significant/insignifant terms each K trials    %
%              if Em.select(i,1) = 0 -> term i insignificant in trial 1   %
%              if Em.select(i,4) = 1 -> term i significant in trial 4     %
%   Xy         Structure: X/Y samples and bootstrap information           %
%    .R        scalar: number of random samples each bootstrap trial      %
%    .X_n      Nxd matrix: normalized parameter vectors                   %
%              X_n(i,1:d) = (X(i,1:d)-X_min)./(X_max-X_min); i = 1,..,N   %
%    .X_min    1xd vector: min values each X column ( = input )           %
%    .X_max    1xd vector: max values each X column ( = input )           %
%    .y        Nx1 vector: Y values supplied by user                      %
%    .id       RxK matrix: index R samples of X for each bootstrap trial  %
%   RT         1xK vector: CPU time (in sec) to construct each emulator   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  IMPORTANT NOTICE                                                       %
%   This code differs from HDMR in that B2 and B3 are not computed ahead  %
%   of time during the initialization phase, but rather computed each     %
%   time their values are needed in 2nd order, 3rd order and refitting.   %
%   This avoids OUT OF MEMORY problems in HDMR with B3 when d > 25.       %
%   Furthermore, the option of backfitting for the 2nd and 3rd order is   %
%   disabled in HDMR_H to speed-up as much as possible emulator           %
%   construction and thus sensitivity analysis for d > 25.                %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  BUILT-IN CASE STUDIES                                                  %
%   Example 1   Multivariate normal benchmark study                       %
%   Example 2   Multivariate normal benchmark study with Σ ≠ identity     %
%   Example 3   Multivariate normal benchmark study with Σ ≠ identity     %
%   Example 4   Ishigami function                                         %
%   Example 5   Ishigami function with Σ ≠ identity                       %
%   Example 6   Function with Σ ≠ identity                                %
%   Example 7   Sobol g function                                          %
%   Example 8   Multivariate normal with Σ ≠ identity                     %
%   Example 9   Example 1 of Chastaing et al. (2012)                      %
%   Example 10  Example 2 of Chastaing et al. (2012)                      %
%   Example 21  Soil temperature modeling                                 %
%   Example 22  Rainfall runoff modeling: hmodel                          %
%   Example 23  Rainfall runoff modeling: SAC-SMA model                   %
%   Example 24  Dynamic foodweb: One-predator-one-prey model              %
%   Example 25  Dynamic foodweb: Two-predators-two-preys model            %
%   Example 26  Dynamic foodweb: Two-predators-two-preys model real data  %
%   Example 27  Simple multivariate function                              %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  MAIN REFERENCE                                                         %
%   Li, G. H. Rabitz, P.E. Yelvington, O.O. Oluwole, F. Bacon,            %
%       C.E. Kolb, and J. Schoendorf (2010), Global sensitivity analysis  %
%       for systems with independent and/or correlated inputs, Journal of %
%       Physical Chemistry A, Vol. 114 (19), pp. 6022 - 6032, 2010        %
%   Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic Sensitivity   %
%       Analysis With Dependent Variables: Covariance-Based Decomposition %
%       of Hydrologic Models, Water Resources Research, 59 (4),           %
%       e2022WR0328346, https://doi.org/10.1029/2022WR032834              %
%                                                                         %
%  MATLAB CODE                                                            %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
%  Version 1.0    January 2018                                            %
%  Version 2.0    July 2024                                               %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% First setup HDMR_H toolbox (add to path)
addpath(pwd,[pwd,'\miscellaneous']); 
% Then go to an example directory, say example_4
cd example_4
% Now one can execute this example by typing: example_4
