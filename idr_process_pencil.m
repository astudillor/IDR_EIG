function [W,C,H,oms] = idr_process_pencil(A,s,m, w0,P,oms_i);
% Computes:
% IDR Factorization (Generalized Hessenberg Relation)
% with orthogonalization
% AW_{m}C_m = W_{m+1}H
%

n = size(A,1);
if nargin < 3
    disp('not enough parameters')
    W = []; H = []; oms = [];
    return
end
if nargin < 4 || isempty(w0)
    w0 = rand(n,1);
    w0 = w0/norm(w0);
end
if nargin < 5 || isempty(P)
    P = rand(n,s); P=P';
end

if nargin < 6 || isempty(oms_i) 
    j = floor(m/(s+1));
    oms_i = ((trace(A)/n) + .0001*randn(j,1) );
end
oms = oms_i;

W = zeros(n,m+1);
H = zeros(m+1,m);
C = eye( m );
M = zeros(s);
f = zeros(s,1);
beta = zeros(s,1);
w = zeros(n,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Create G_0
W(:,1) = w0;
for i=1:s,
    w = A*W(:,i);
    for j=1:i
        H(j,i) = w'*W(:,j); 
        w = w - H(j,i)*W(:,j);
    end
    H(i+1,i) = norm(w);
    W(:,i+1) = w/H(i+1,i);
    M(:,i) = P*W(:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k = s+1:m,
    if mod(k,s+1) == 0,
        %Enter to the G_j
        om = oms(ceil((k+1)/(s+1)) -1);
        i = 0;
    end
    w = W(:,k);
    f = (P*w);
    beta = M\f;
    w = w - W(:,(k-s):(k-1))*beta;
    w =  A*w - om*w; %(A-om I)v

    H((k-s):k-1,k) = -om*beta;
    H(k,k) = om; 
    H(k+1,k) = 1; 

    %DGKS orthogonalization
    %for t=1:2    
        c = W(:,k-i+1:k)'*w;
        w = w - W(:,k-i+1:k)*c;
        H(k-i+1:k,k) = H(k-i+1:k,k) + c;
    %end

    H(k+1,k) = norm(w); 
    w = w/H(k+1,k);

    C((k-s):k-1,k) = -beta;
    W(:,k+1) = w; 
    M(:,1:s-1) = M(:,2:s);
    M(:,s) =  f;
    i = i + 1;   %Another vector in G_j
end