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
        
    case 'ordered'
        opt                     = genetic.tools.keyValuePairs(varargin{:});
        [t,X,iCanon,canonNames] = genetic.plot.nDimOrdered(n, opt.N);
        % To close the disk, the first element are repeated
        t                       = [t, t(1)];
        X                       = [X, X(:,1)];
        %
        r                       = linspace(0,opt.rmax, opt.nRings);
        Y                       = zeros(length(t),opt.nRings);
        x1                      = Y;
        x2                      = Y;        
        for i = 1:length(t)
            xi = X(:,i);
            for j = 1:opt.nRings
                Y(i,j) = fun(x0 + r(j) * xi);
                x1(i,j) = r(j) * cos(t(i));
                x2(i,j) = r(j) * sin(t(i));
            end
        end
        surf(x1,x2, Y);
        % Adding canonical directions
        hold on
        for i = 1:length(iCanon)
            plot3(x1(iCanon(i),:),x2(iCanon(i),:),Y(iCanon(i),:),'r','lineWidth',2)
            h = text(1.3*x1(iCanon(i),end), 1.3*x2(iCanon(i),end),Y(iCanon(i),end),canonNames{i});
            set(h,'color','red');
        end
end
end