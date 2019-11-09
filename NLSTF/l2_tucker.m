function  A   =  l2_tucker( D, X,lambda )

A        =    zeros( size(D{1},2),size(D{2},2),size(D{3},2),size(X,4));




DTX = TensorChainProductT(X,D,1:3);

wt=1;
N = numel(D);
for i = ndims(X) : -1 : N+1
    wt = ones(size(X,i),1)*wt';
    wt = wt(:);
end

mu        =  lambda;

for i = 3 : -1 : 1
    [U,S] = eig(D{i}'*D{i});
    P{i} = U;
    wt = diag(S)*wt';
    wt = wt(:);
end
wt = reshape(wt+mu,size(A));  % diagonal vector of (ita*kron(S_1,..S_n) + 2*gamma*I)^(-1)
Z = TensorChainProductT(DTX,P,1:3);
Z = Z./wt;
A = TensorChainProduct(Z,P,1:3);




