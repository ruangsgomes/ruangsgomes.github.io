% Clear memory and close matlab Figures
clc; clear all; close all;

%ndega =     3
% EPR_out_MonteCarlo_MODEGA(:,1).bestpop
% ans =
%         0         0    2.0000         0         0         0         0
%         0         0         0    1.0000    1.0000    1.0000         0
%    1.0000    0.2000         0         0         0         0    1.0000

% vector of exponents
expo = [-2:0.5:2];
% optimum number of polynomial terms
ndega = 3;
% best EPR parameters
best_par = rand(ndega+1,1);
best_par(end) = 14e-8;
% input variables
X = rand(200,5);
% "random" best population
pop = zeros(ndega,size(X,2));
for i = 1:ndega
    pop(i,:) = datasample(expo,size(X,2));
end


% print final equation
print_EPR_equation (X,best_par,pop);

