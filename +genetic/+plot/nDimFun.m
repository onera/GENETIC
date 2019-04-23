% GENETIC.PLOT.NDIMFUN - Plot scalar function in some neighbourhood in Rn
%
% Syntax
%   genetic.plot.nDimFun(x0, fun, param, varargin)
%
% Parameters
%   x0      : point in Rn around which the space is sampled
%   fun     : objective function 
%   param   : name of the sampling of Rn to be used. Possible values are:
%             'polar', 'ordered'
%   varargin : key/value pairs for the options,
%              - 'N': number of sampling points per ring 
%              - 'nRings': number of rings,starting from x0
%              - 'rmax': maximum radius starting from x0 
%              - 'cst': constraints structure to that infeasible points are
%              shown
%
% References
%   W.L. Goffe, Visualizing Multi-Dimensional Functions in Economics. 1999

% Copyright 2018 ONERA
%
% This file is part of the GENETIC project.
%
% GENETIC is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License version 3 as
% published by the Free Software Foundation.
%
% GENETIC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with GENETIC.  If not, see <https://www.gnu.org/licenses/lgpl-3.0>.
function nDimFun(x0, fun, param, varargin)
n           = length(x0);
if ~isempty(varargin)
    opt     = genetic.tools.keyValuePairs(varargin{:});
else
    opt     = struct();
end
default_opt = struct('N',30,'nRings',2,'rmax',1,'cst',[]);
opt         = genetic.tools.fillEmptyFields(opt, default_opt);
% constraints
cst         = genetic.constraints(opt.cst);

switch param
    case 'polar'
        % See the article:
        %   W.L. Goffe, Visualizing Multi-Dimensional Functions in
        %   Economics. 1999
        [t,X]   = genetic.plot.nDimPolar(n, opt.N);
        r       = linspace(opt.rmax/opt.nRings, opt.rmax, opt.nRings);
        %
        Y = zeros(opt.N,opt.nRings);
        x1 = Y;
        x2 = Y;
                feasible                = Y == 0;

        for i = 1:opt.N
            xi = X(:,i);
            for j = 1:opt.nRings
                xij     = x0 + r(j) * xi;
                Y(i,j)  = fun(xij);
                x1(i,j) = r(j) * cos(t(i));
                x2(i,j) = r(j) * sin(t(i));
                feasible(i,j)   = cst.satisfied(xij);

            end
        end
                % Complete the data with the value at x0
        x1  = [0*x1(:,1), x1];
        x2  = [0*x2(:,1),x2];
        y0  = fun(x0);
        
        Y   = [y0*ones(length(t),1),Y];
        feasible = [cst.satisfied(x0) * ones(length(t),1)>0,feasible];
        %
        figure
        surf(x1,x2, Y);
        hold on
        % Plotting infeasible points
        plot3(x1(~feasible),x2(~feasible),Y(~feasible),'m.','MarkerSize',20);

    case 'ordered'
        [t,X,iCanon,canonNames] = genetic.plot.nDimOrdered(n, opt.N);
        % To close the disk, the first element are repeated
        t                       = [t, t(1)];
        X                       = [X, X(:,1)];
        %
        r                       = linspace(opt.rmax/opt.nRings, opt.rmax, opt.nRings);
        Y                       = zeros(length(t),opt.nRings);
        x1                      = Y;
        x2                      = Y;
        feasible                = Y == 0;
        for i = 1:length(t)
            xi = X(:,i);
            for j = 1:opt.nRings
                xij             = x0 + r(j) * xi;
                Y(i,j)          = fun(xij);
                x1(i,j)         = r(j) * cos(t(i));
                x2(i,j)         = r(j) * sin(t(i));
                feasible(i,j)   = cst.satisfied(xij);
            end
        end
        % Complete the data with the value at x0
        x1  = [0*x1(:,1), x1];
        x2  = [0*x2(:,1),x2];
        y0  = fun(x0);
        
        Y   = [y0*ones(length(t),1),Y];
        feasible = [cst.satisfied(x0) * ones(length(t),1)>0,feasible];
        %
        figure
        surf(x1,x2, Y);
        hold on
        % Adding canonical directions
        for i = 1:length(iCanon)
            plot3(x1(iCanon(i),:),x2(iCanon(i),:),Y(iCanon(i),:),'r','lineWidth',2)
            h = text(1.3*x1(iCanon(i),end), 1.3*x2(iCanon(i),end),Y(iCanon(i),end),canonNames{i});
            set(h,'color','red');
        end
        % Plotting infeasible points
        plot3(x1(~feasible),x2(~feasible),Y(~feasible),'m.','MarkerSize',20);
        %
%         plotInfeasibility(x1,x2,Y,feasible);
end
end

% function plotInfeasibility(x1, x2, Y, feasible)
% for i = 1:size(feasible,1)-1
%     for j = 1:size(feasible,2)-1
%         if ~feasible(i,j)
%             ip1 = i + 1;
%             for j2 = j-1:j+1
%                 if ~feasible(ip1,j2)
%                     plot3([x1(i,j),x1(ip1,j2)],[x2(i,j),x2(ip1,j2)],[Y(i,j),Y(ip1,j2)],'k');
%                 end
%             end
%         end
%     end
% end
% end