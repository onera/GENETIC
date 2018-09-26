function [t,X,iCanon,canonName] = nDimOrdered(n, N, order)
if nargin < 3
order = 1:n;
end
r       = pi/n;
%
tCanon  = zeros(2*n+1,1);
canonName = {};
for i = 1:2*n+1
    tCanon(i) = (i-1) * r;
    if i < 2*n+1
        if i <= n 
            idx = sprintf('%d',order(i));
            pre = '';
        elseif i < 2*n
            idx = sprintf('%d',order(rem(i,n)));
            pre = '-';
        else
            idx = sprintf('%d',order(n));
            pre = '-';
        end
        canonName{i} = [pre,'x_',idx];
    end
end
Nleft       = N - 2*n;
nPerQuadran = ceil(Nleft / (2*n));
t           = [];
iCanon      = zeros(2*n,1);
for i = 1:2*n
    iCanon(i)   = max([length(t),1]);
    t           = [t(1:end-1), linspace(tCanon(i), tCanon(i+1),nPerQuadran + 2 )];
end
t(end) = [];
%
N = length(t);
X = zeros(n,N);
I = eye(n);


for i = 1:N
    if t(i) == 0 
        quadran = 1;
    else
        quadran = ceil(t(i) / r);
    end
    a1 = 1;
    a2 = 1;
    if quadran < n 
        id1 = quadran;
        id2 = id1 + 1;
    elseif quadran == n 
        id1 = quadran;
        id2 = 1;
        a2  = -1;
    elseif quadran == 2 * n
        id1 = n;
        id2 = 1;
        a1  = -1;
    else
        id1 = rem(quadran, n);
        id2 = id1 + 1;
        a1 = -1;
        a2 = -1;
    end
    alpha   = cos((t(i) - (quadran-1) * r) * n / 2);
    e1      = a1 * I(:,order(id1));
    e2      = a2 * I(:,order(id2));
    X(:,i)  = alpha * e1 + (1 - alpha) * e2;
    X(:,i)  = X(:,i) / norm(X(:,i));
end



end