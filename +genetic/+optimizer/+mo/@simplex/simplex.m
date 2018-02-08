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
classdef simplex < genetic.optimizer.mono & genetic.optimizer.simpleScheme
   % SIMPLEX - Nelder-Mead simplex
   %
   %  The Nelder-Mead simplex [1] is a derivative-free method
   %  relying on the evolution of a simplex in the search-space.
   %
   %  At each iteration, one vertex of the simplex is moved according
   %  to some mechanism: reflection, expansion or contraction. Those 
   %  mechanisms enables to guide the research towards interesting
   %  directions.
   %
   %  Note that it has been proven that the Nelder-Mead Simplex can fail to
   %  converge even on convex functions but it has also proven to be fairly
   %  efficient on other cases. Besides, restarting the simplex is
   %  generally enough to avoid beeing stuck in a non-stationary point.
   %  Note that te Nerlder-Mead simplex is not a generating set search
   %  method.
   %
   %  Parameters
   %     tolSimplexSize     : tolerance on the size of the simplex below which
   %                          the algorithm stops (1e-8)
   %     restartSimplex     : number of allowed restart for the simplex (0)
   %     initialSimplex     : nature of the initial simplex. Can be 'regular',
   %                          'right' ('regular')
   %     initialSimplexSize : size of the initial simplex (default is
   %                          determined depending on x0)
   %     reflectCoef        : coefficient for the reflection phase (1)
   %     expandCoef         : coefficient for the expand phase (2)
   %     contractCoef       : coefficient for the contraction phase (1/2)
   %     
   % References
   %     [1] J.A. Nelder and R. Mead, A Simplex Method for Function Minimization,
   %         The computer Journal, 1965. 
   %
   properties
      restartSimplex       = 0;
      tolSimplexSize       = 1e-8;
      %
      initialSimplex       = 'regular';
      initialSimplexSize   = [];
      %
      reflectCoef          = 1;
      expandCoef           = 2;
      contractCoef         = 1/2;
   end
   properties (SetAccess=protected)
      simplexX
      simplexY
      simplexCentroid
      phase           = 'reflect';
      reflectedX
      reflectedY
   end
   methods
      %% constructor
      function self = simplex(varargin)
         self@genetic.optimizer.mono(varargin{:});
         self@genetic.optimizer.simpleScheme();
         %
         self.popSize      = 1;
         %
         self.methodName   = 'simplex';
         self.longName     = 'Nelder-Mead Simplex';
         self.toDisplay{2}.name = 'iter';
         % Adding stopping criterion concerning the size of the simplex
         simplexTooSmall = struct('name', 'simplexTooSmall', 'cleanStop', true, 'condByConstraints', false);
         self.addStopTest(simplexTooSmall);
         % Adding printing function
         simpSize        = struct('name','simplexSize','call','printSimplexSize','align','c','dim',16);
         self.addDisplayElement(simpSize);
      end
      %% updateBarycenter
      function updateBarycenter(self)
         self.simplexCentroid = 1/self.xDim * sum(self.simplexX(:,1:end-1),2);
      end
      %% sortSimplex
      function sortSimplex(self)
         [self.simplexY, idx] = sort(self.simplexY,'ascend');
         self.simplexX        = self.simplexX(:,idx);
      end
      %% simplexSize
      function s = simplexSize(self)
         v1 = self.simplexX(:,1);         
         s  = norm(self.simplexX(:,2:self.xDim + 1) - repmat(v1,1,self.xDim),1) / max(1, norm(v1,1));
      end
      %% evolve1
      function group = evolve(self, group)
         switch self.phase
            case 'reflect'
               % Computation of reflected point
               newX              = self.simplexCentroid + self.reflectCoef * (self.simplexCentroid - self.simplexX(:,end));
               self.reflectedX   = newX;
            case 'expand'
               newX              = self.simplexCentroid + self.expandCoef*(self.reflectedX -self.simplexCentroid);
            case 'contract'
               newX              = self.simplexX(:,end) + self.contractCoef * (self.simplexX(:,end) - self.simplexCentroid);
            case 'shrink'
               x1                = self.simplexX(:,1);
               for i = 2:self.xDim + 1
                  self.simplexX(:,i) = x1 + 1/2 * (self.simplexX(:,i) - x1);
               end
               group.spawnIndividual(self.xDim-1);
               newX              = self.simplexX(:,2:end);
            case 'restart'
               x1                = self.simplexX(:,1);
               self.simplexX     = genetic.optimizer.mo.simplex.createSimplex('regular', x1);
               group.spawnIndividual(self.xDim-1);
               newX              = self.simplexX(:,2:end);
               %
               self.restartSimplex = self.restartSimplex -1;
         end
         %
         group.moveTo(newX)
      end

      %% 
      function postEval(self, group)
         if self.nGen > 0
            switch self.phase % the transformation that has just been made
               case 'reflect'
                  self.reflectedY   = group.getObjective();
                  fxn               = self.simplexY(:,end-1);
                  fx1               = self.simplexY(:,1);
                  if fx1 <= self.reflectedY && self.reflectedY < fxn
                     self.replaceInSimplex(self.xDim + 1, self.reflectedX, self.reflectedY);
                     nextPhase   = 'reflect';
                  elseif self.reflectedY < fx1
                     nextPhase   = 'expand';
                  else
                     nextPhase   = 'contract';
                  end
               case 'expand'
                  fe = group.getObjective();
                  if fe < self.reflectedY
                     self.replaceInSimplex(self.xDim + 1, group.getValue(), fe);
                  else
                     self.replaceInSimplex(self.xDim + 1, self.reflectedX, self.reflectedY);
                  end
                  nextPhase = 'reflect';
               case 'contract'
                  fc = group.getObjective();
                  if fc < self.simplexY(end)
                     self.replaceInSimplex(self.xDim + 1, group.getValue(), fc);
                     nextPhase = 'reflect';
                  else
                     nextPhase = 'shrink';
                  end
               case 'shrink'
                  self.simplexY(2:end) = group.getObjective();
                  [~,id]               = min(group.getObjective());
                  group.killIndividuals(setdiff(1:group.len,id));
                  nextPhase = 'reflect';
               case 'restart'
                  self.simplexY(2:end) = group.getObjective();
                  while group.len > 1
                     group.at(2).kill();
                  end
                  nextPhase = 'reflect';
            end
            self.phase = nextPhase;
         else
            self.simplexX        = group.getValue();
            [self.simplexY, idx] = sort(group.getObjective);
            self.simplexX        = self.simplexX(:,idx);
            % kill all individuals except 1
            while group.len > 1
               group.at(2).kill();
            end
         end
         % Resort simplex values
         self.sortSimplex();
         % Update simplex barycenter
         self.updateBarycenter();
         % Compute simplex size
         if self.simplexSize < self.tolSimplexSize
            self.phase = 'restart';
         end
      end
      %%
      function replaceInSimplex(self, i, x, y)
         self.simplexX(:,i)   = x;
         self.simplexY(i)     = y;
      end
      %%
      function group = makeGroup(self)
         group = genetic.population.group(self.xDim + 1);
      end
      %% initPopulation
      function x = initPopulation(self, x0)
         %
         if isempty(x0)
            x0 = genetic.population.group.initValue(self.initMethod, 1, [], self.xMin, self.xMax, self.regulariseInit);
         end
         %
         if isa(self.initialSimplex,'double')
            x = self.initialSimplex;
            return
         end
         %
         scale = [];
         if ~isempty(self.initialSimplexSize)
            % If the initial size is provided
            scale = self.initialSimplexSize;
         else
            % Otherwise, try to determine some size according to the
            % bounding box (if any)
            if self.isBoxBounded
               scale = min(norm(self.xMax - x0,inf),norm(self.xMin - x0,inf));
            end
            
         end
         %
         x = genetic.optimizer.mo.simplex.createSimplex(self.initialSimplex, x0, scale);
      end
      %% simplexTooSmall
      function [out, msg] = simplexTooSmall(self, group)
         s     = self.simplexSize;
         out   = s < self.tolSimplexSize & self.restartSimplex == 0;
         msg   = '';
         if out 
            msg = sprintf('Size of the simplex below tolerance: %1.2e (tol: %1.2e)', s, self.tolSimplexSize);
         end
      end
      %%
      function out = printSimplexSize(self, group)
         out = sprintf('%1.5e',self.simplexSize);
      end
   end
   %%
   methods (Static)
      function v = createSimplex(key, x0, scale)
         % Creating an initial simplex of given size
         if nargin < 3 || isempty(scale)
            scale = max(norm(x0,inf),1);
         end
         x0       = x0(:);
         n        = length(x0);
         %
         v        = [x0,eye(n,n)];
         switch key
            case 'regular'
               % Regular simplex - all edges have same length.
               p           = sqrt(n + 1) + n -1;
               q           = sqrt(n + 1) - 1;
               alpha       = scale/(n * sqrt(2)) * [p, q];
               v(:,2:n+1)  = (x0 + alpha(2) * ones(n,1)) * ones(1,n);
               for i = 2:n+1
                  v(i-1,i) = x0(i-1) + alpha(1);
               end
            case 'right'
               % Right-angled simplex
               alpha       = scale*ones(n+1,1);
               for i = 2:n+1
                  v(:,i)   = x0 + alpha(i) * v(:,i);
               end
         end
      end
      %
      function plotSimplex(x,y)
         n = length(y);
         for i = 1:n
            cur = x(:,i);
            if i < n
               next = x(:,i+1);
            else
               next = x(:,1);
            end
            plot([cur(1), next(1)],[cur(2),next(2)],'b.-')
            hold on
            
         end
         [~,id] = min(y);
         plot(x(1,id),x(2,id),'ro')
         
      end
   end
end