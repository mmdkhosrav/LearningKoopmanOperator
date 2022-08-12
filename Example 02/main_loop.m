%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% This code attempts to solve the following learning problem  
%       min_K   E(K) + lambda * R(K),
%       s.t.    K in C,
% where E(.) is the sum squared error loss function, C is a given constraint, 
% R(.) is the regularization term. For C and R, we have the following cases:
%   1. R: square of operator norm of K, i.e.,  R(K) = ||K||^2 
%   2. R: square of Frobenius norm of K, i.e.,  R(K) = ||K||_F^2 
%   3. R: nuclear norm of K, i.e.,  R(K) = ||K||_* 
%   4. C: rank of K, i.e.,  C = {K | rank(K) <= r} 
%   5. R&C: if R=0 and C = L(K) (the whole space of bounded operators), the
%      program is equivalent to the EDMD approach
% For generating the data, two trajectories of the Van der Pol system is 
% considered. More details are provided in Example 2 in the reference above.                                                   
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
close all;  clear all;  clc;
%--------------------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('Learnin Koopman Operator:                                          ')
disp('Ref: Representer Theorem for Learning Koopman Operators            ') 
disp('Link: https://arxiv.org/abs/2208.01681                             ')
disp('                                                                   ')
disp('This code attempts to solve the following learning problem         ') 
disp('       min_K   E(K) + lambda * R(K),                               ')
disp('       s.t.    K in C,                                             ')
disp('where E(.) is the sum squared error loss function, C corresponds to')
disp('given constraint, R(.) is the regularization term. For C and R, we ')
disp('consider the following cases:                                      ')
disp('   1. R: square of operator norm of K, i.e., R(K) = ||K||^2        ') 
disp('   2. R: square of Frobenius norm of K, i.e., R(K) = ||K||_F^2     ') 
disp('   3. R: nuclear norm of K, i.e.,  R(K) = ||K||_*                  ') 
disp('   4. C: rank of K, i.e.,  C = {K | rank(K) <= r}                  ') 
disp('   5. R&C: if R = 0 and C = L(K) (the space of bounded operators), ') 
disp('      the program is equivalent to the EDMD approach.              ')
disp('For generating the data, two trajectories of the Van der Pol system')
disp('is considered. More details and discussion are provided in Example ')
disp('2 in the reference above.                                          ')
disp('                                                                   ')
disp('Mohammad Khosravi                                                  ')
disp('Email: mohammad.khosravi@tudelft.nl                                ')
disp('Delft Center for Systems and Control (DCSC)                        ')
disp('Delft University of Technology (TU Delft)                          ') 
disp('August 2022                                                        ')
disp('-------------------------------------------------------------------')
disp('The program is started!')
format shortg
c = clock;
disp(['Start time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')
%--------------------------------------------------------------------------

LW = 'linewidth';       FS = 'FontSize';        MS = 'MarkerSize';
LOC = 'Location';       INT = 'Interpreter';    LX = 'Latex';   
DN = 'DisplayName';     
OLS = 'OutlierSize';    FG = 'FactorGap';

CL = 'Color';
blue_CL     = [0 0.45 0.74];    sky_CL      = [0 0.75 0.75];
olive_CL    = [0.75 0.75 0];    orange_CL   = [0.85 0.33 0.1];
purple_CL   = [0.75 0 0.75];

%--------------------------------------------------------------------------
% adding the path of support files
currentFolder = pwd;
addpath([currentFolder(1:end-10),'support_codes'])

%--------------------------------------------------------------------------
% Monte Carlo loops

% SNR level 10 dB
dBvalx = 10;
parfor MC_num = 1:120
    LKO_comp_fun(dBvalx,MC_num);    
end

% SNR level 20 dB
dBvalx = 20;
parfor MC_num = 1:120
    LKO_comp_fun(dBvalx,MC_num);    
end

% SNR level 30 dB
dBvalx = 30;
parfor MC_num = 1:120
    LKO_comp_fun(dBvalx,MC_num);    
end
















