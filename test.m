clc
clear all
close all
format long

n = 3000;
A = sprandn(n,n,0.2);
eA = sort(eig(full(A)));
s = 8;
m = 100;

w0 = rand(n,1); w0 = w0/norm(w0);

P = rand(s,n);
j = floor(m/(s+1));
oms_i = ((trace(A)/n) + .0001*randn(j,1) );
    
tic
[W1,C1,H1,oms1] = idr_process_pencil(A,s,m, w0,P,oms_i);
eIdr_p = sort(eig(H1(1:m,1:m),C1));
t_idr_p = toc

tic
[W2,H2,oms1] = idr_process_standard(A,s,m, w0,P,oms_i);
eIdr_s = sort(eig(H2(1:m,1:m)));
t_idr_s = toc

tic
[V,H3] = arnoldi(A,w0,m);
eAr = sort(eig(H3(1:m,1:m)));
t_ar = toc


fprintf('|AWC - WH| = %g\n',norm(A*W1(:,1:end-1)*C1 - W1*H1))
fprintf('|AW - WH| = %g\n',norm(A*W2(:,1:end-1) - W2*H2))
fprintf('|AV - VH| = %g\n',norm(A*V(:,1:end-1) - V*H3))

plot(real(eA),imag(eA),'ok')
hold on
plot(real(eAr),imag(eAr),'*b')
plot(real(eIdr_s),imag(eIdr_s),'*r')
plot(real(eIdr_p),imag(eIdr_p),'*g')
legend('Eigenvalues of A', 'Ritz Arnoldi','Ritz IDR standard','Ritz IDR pencil','location','best')