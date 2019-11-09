

function  A   =  Nonnegative_SC( D, X, par )

A        =    zeros( size(D,2), size(X, 2) );
V        =    zeros( size(D,2), size(X, 2) );
T        =    50;
DTD       =   D'*D;
DTX       =   D'*X;
Ek        =   eye(par.K);
mu        =   0.01;
ro        =   1.1;
for  i  =  1:T
    S         =   (DTD + mu*Ek)\(DTX + mu*(A-V/(2*mu)) );
    A         =   max( soft(S+V/(2*mu), par.lambda/(2*mu)), 0);
    V         =   V + 2*mu*( S - A );
    mu        =   mu*ro;
end
