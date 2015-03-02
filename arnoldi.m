function [W,H] = arnoldi(A,w0,s)

n = size(A,1);
w0 = w0/norm(w0);

W = zeros(n,(s+1));
H = zeros(s+1,s);


W(:,1) = w0;
for i=1:s,
    w = A*W(:,i);
    for j=1:i
        H(j,i) = w'*W(:,j); 
        w = w - H(j,i)*W(:,j);
    end
    H(i+1,i) = norm(w);
    W(:,i+1) = w/H(i+1,i);
end