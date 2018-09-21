function nDimFun(x0, fun, param, varargin)
n = length(x0);
switch param
    case 'polar'
        % See the article:
        %   W.L. Goffe, Visualizing Multi-Dimensional Functions in
        %   Economics. 1999
        opt     = genetic.tools.keyValuePairs(varargin{:});
        [t,X]   = genetic.plot.nDimPolar(n, opt.N);
        r       = linspace(0,opt.rmax, opt.nRings);
        %
        Y = zeros(opt.N,opt.nRings);
        x1 = Y;
        x2 = Y;
        for i = 1:opt.N
            xi = X(:,i);
            for j = 1:opt.nRings
                Y(i,j) = fun(x0 + r(j) * xi);
                x1(i,j) = r(j) * cos(t(i));
                x2(i,j) = r(j) * sin(t(i));
            end
        end
        surf(x1,x2, Y);
        
end
end