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
%
classdef annealing < genetic.optimizer.mono & genetic.optimizer.simpleScheme
   % ANNEALING - Simulated annealing
   %
   % Description
   %
   % Parameters
   %  T0: initial temperature (estimated by default)
   %  Tmin: minimal temperature (1e-6)
   %  collingScheme: scheme for the cooling procedure
   %                 'exp' (default), 'log', 'linear', 'fast'
   %  maxIterAtT: maximum number of iteration at a temperature (30)
   %
   % References
   % [1] S. Kirkpatrick, C.D. Gelatt and M.P. Vecchi, Optimization by 
   %     Simulated Annealing, Science, 1983.
   %
   properties
      % Initial temperature
      nEvalsT0       = 50;
      probaStart     = 0.7;
      T0
      % Cooling scheme
      Tmin           = 1e-6;
      coolingScheme  = 'exp';
      maxIterAtT     = 100;
      alpha          = 0.95;
      %
      neigh          = 'default';
      delta          = 0.1;
      %
   end
   properties(SetAccess=protected)
      annealingParams = [];
      % State
      T
      nInner         = 0;
      nOuter         = 0;
      x 
      y
      %
      computeT0      = false;
      x0             = [];
      %
      nAccept        = 0;
      
   end
   methods
      %% constructor
      function self = annealing(varargin)
         self@genetic.optimizer.mono(varargin{:});
         self@genetic.optimizer.simpleScheme();
         self.methodName         = 'annealing';
         self.longName           = 'Simulated annealing';
         self.needFiniteBounds   = true;
         self.popSize            = 1;
         self.toDisplay{2}.name  = 'iter';

%          self.printOnlyImp       = true;
%          self.modPrintIter       = 10;
         %
         self.computeT0          = isempty(self.T0);
         % Adding printing function
         temp        = struct('name','Mean temperature','call','printTemperature','align','c','dim',16);
         self.addDisplayElement(temp);
      end
      %%
      function out = printTemperature(self, dummy)
         out = sprintf('%.5e',mean(self.T));
      end
      %%
      function group = makeGroup(self)
         if ~self.computeT0
            group    = makeGroup@genetic.optimizer.mono(self);
            self.T   = self.T0 * ones(self.xDim, 1);
         else
            group = genetic.population.group(self.nEvalsT0);
         end
      end
      %%
      function x = initPopulation(self, x0)
         x        = initPopulation@genetic.optimizer.mono(self, x0);
         if self.computeT0
            self.x0     = x;
            x           = zeros(self.xDim, self.nEvalsT0);
            x(:,1)      = genetic.population.group.initValue(self.initMethod, 1, [], self.xMin, self.xMax, false, self.defaultInfinity);
            %
            tmp         = self.neigh;
            self.neigh  = 'default';
            for i = 1:self.nEvalsT0 -1
               self.x   = x(:,i);
               x(:,i+1) = self.neighbourhood();
            end
            self.neigh = tmp;
         end
      end
      %%
      function postEval(self, group)
         if self.nGen  == 0 && self.computeT0
            y0    = group.getObjective();
            fSum  = 0;
            nImp  = 0;
            for i = 2:length(y0)
               if y0(i-1) < y0(i)
                  nImp = nImp + 1;
                  fSum = fSum + (y0(i) - y0(i-1));
               end
            end
            if nImp >0
               fSum = fSum/nImp;
            end
            self.T0     = -fSum ./ (log(self.probaStart));
%             self.T0     = abs(fSum / log(1/self.probaStart - 1));
            self.T0     = self.T0 * ones(self.xDim,1);
            self.T      = self.T0;
            group.killIndividuals(2:group.len);
            group.moveTo(self.x0)
         end
      end
      %% evolve
      function group = evolve(self, group)
         self.nInner = self.nInner + 1;
         if self.nInner == 1 && all(self.T == self.T0)
            self.x    = group.getValue();
            self.y    = group.getObjective();
         else
            % New candidate solution
            candidatex     = group.getValue();
            candidatey     = group.getObjective();
            % Acceptation
            delta_energy   = (candidatey - self.y)';
            downHill       = delta_energy < 0;
            upHill         = ~downHill;
            accept         = downHill | (upHill & (exp(-delta_energy ./ self.T)> rand(group.len,1)));
            % Replacement of current solution
            if accept
               self.nAccept   = self.nAccept + 1;
               self.x         = candidatex;
               self.y         = candidatey;
            end
%             newX = self.x;
%             newX(:,accept) = candidatex(:,accept);
%             self.x         = newX;
%             newY           = self.y;
%             newY(:,accept) = candidatey(:,accept);
%             self.y         = newY;
         end
         
%          if self.nAccept == self.reAnnealing
%             s                          = 
%             self.annealingParameters   = log( self.T0 ./ self.T .* max(s) ./s);
%          end
         
         if self.nInner >= self.maxIterAtT
            % Decreasing temperature
            self.nInner = 0;
            self.nOuter = self.nOuter + 1;
            self.decreaseTemperature();
         end
         % Moving to a neighbourhood of the actual point
         group.moveTo(self.neighbourhood(group));
      end
      %% 
      function T = decreaseTemperature(self)
         if isa(self.coolingScheme,'function_handle')
            self.T = self.coolingScheme(self.nOuter, self.T);
            return
         end
         switch self.coolingScheme
            case 'exp'
               T = self.alpha * self.T;
            case 'log'
               T = self.T0/log(self.nOuter);
            case 'fast'
               T = self.T0/(self.nOuter);
         end
         self.T = T;
      end
      %%
      function xn = neighbourhood(self, g) 
         if isa(self.neigh,'function_handle')
            xn = self.neigh(self.x, self.T, self.xMin, self.xMax);
            xn = self.constraints.projectOnBounds(xn, self.xMin, self.xMax);
            return
         end
         D = self.xMax - self.xMin;
         switch self.neigh
            case 'default'
               xn = self.x + self.delta * D .* randn(self.xDim, self.popSize);
            case 'linear'
               xn = self.x + self.delta * max(self.T,1) .* randn(self.xDim, self.popSize);
            case 'perm'
               xn = self.x + self.delta * (randperm(self.xDim)'==self.xDim)*randn(1);
            case 'asa'
               u  = randn(self.xDim,1);
               yi = sign(u - 1/2) .* self.T .* ((1+1./self.T).^(abs(2*u -1)) -1 );
               xn = self.x + yi .* D;
            otherwise
               error('Unknown neighbourhood type')
         end
         xn = self.constraints.projectOnBounds(xn, self.xMin, self.xMax);
      end

      function [out, msg] = stop(self, group)
         out = self.T <= self.Tmin;
         msg = '';
         if out 
            msg = sprintf('Minimum temperature reached: %1.2e (%1.2)',self.T, self.Tmin);
         end
      end
   end
end