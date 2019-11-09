function  A   =  sparse_tucker( D, X, par )

A        =    zeros( size(D{1},2),size(D{2},2),size(D{3},2),size(X,4));
V        =    zeros( size(A) );
T        =    50;


DTX = TensorChainProductT(X,D,1:3);

wt=1;
N = numel(D);
for i = ndims(X) : -1 : N+1
    wt = ones(size(X,i),1)*wt';
    wt = wt(:);
end

mu        =   0.01;

for i = 3 : -1 : 1
    [U,S] = eig(D{i}'*D{i});
    P{i} = U;
    wt = diag(S)*wt';
    wt = wt(:);
end
wt = reshape(wt+mu,size(A));  % diagonal vector of (ita*kron(S_1,..S_n) + 2*gamma*I)^(-1)
Z = TensorChainProductT(DTX,P,1:3);
Z = Z./wt;
bbb = TensorChainProduct(Z,P,1:3);



% aaa= (DTD + mu*Ek)\eye(size(D,2));
% bbb=aaa*DTX;


for  i  =  1:T
    Z = TensorChainProductT(mu*(A-V/(2*mu)),P,1:3);
Z = Z./wt;
ccc = TensorChainProduct(Z,P,1:3);
    S         =   bbb+ccc ;  
 A         =   soft(S+V/(2*mu), par/(2*mu));
    V         =   V + 2*mu*( S - A );


%     fun(i)    =   0.5*sum(sum((X-D*A).^2)) + par.lambda*sum(sum(abs(A)));
end
