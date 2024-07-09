clear all
close all
clc


A=  [1.0, 4.0, 6.0, 9.0,
  -3.0, 2.0, 7.0, 7.0, 
  2.0, 0.0, 3.0, -2.0, 
  5.0, 1.0, -3.0, 4.0 ];

% rng(42);
LIM = 32767;
n=4;
% vpa(50)
while 1
    
    
    scale = LIM * 2;  % 缩放因子
    shift = -LIM;  % 平移因子

    random_matrix = (rand(n, n) * scale) + shift;

%     A = rand(n,n);
%     A = random_matrix;

%     save A
%     load A1
    eig(A)

    [H,WR,WI] = dgeev(A);
    
    eig_val = (WR+ 1j.*WI)
    
    if( sum(abs(eig(A)))-sum(abs( eig_val)) > 1e-4)
        error('123');
    end
    if any(isnan(eig_val))
        error('456');
    end
    pause(1);
end
