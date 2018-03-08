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
classdef constraints < handle
  
   properties
      A        = [];
      b        = [];
      Aeq      = [];
      beq      = [];
      c        = [];
      ceq      = [];
      xMin     = -inf;
      xMax     = inf;
      %
      skip     = {};
      %
      tolCst   = 1e-8;
      %
      lastViolation = [];
   end

   
   methods
      %% constraints
      function self = constraints(cstStructure)
         % constructor of the constraints class
         if ~isempty(cstStructure)
            fieldsNames = fieldnames(cstStructure);
            for i = 1:length(fieldsNames)
               try
                  self.(fieldsNames{i}) = cstStructure.(fieldsNames{i});
               catch
                  fprintf('''%s'' is not an attribute of %s, discarding.\n',fieldsNames{i},class(self));
               end
            end
         end
         %
      end
      %% empty
      function out = empty(self)
         % test whether these constraints are empty
         out = isempty(self.A) && isempty(self.Aeq) && isempty(self.c) && isempty(self.ceq) && all(self.xMin==-inf) && all(self.xMax==inf);
      end
      %% penalise
      function [fun, n] = penalise(self, toPenalise)
         % returns a penalty functions f(x, sigma) where sigma is the
         % coefficient of the penalty
         penalisationStr = '';
         n = 0;
         if ~isempty(strcmp('li',toPenalise)) && ~isempty(self.A)
            A                 = self.A;
            b                 = self.b;
            penalisationStr   = [penalisationStr, ' + sigma *max(max(A * x - b),0)'];
            n = n+1;
         end
         %
         if ~isempty(strcmp('le',toPenalise)) && ~isempty(self.Aeq)
            Aeq               = self.Aeq;
            beq               = self.beq;
            penalisationStr   = [penalisationStr,'+ sigma * norm(Aeq * x - beq,2)^2 '];
            n                 = n+1;
         end
         %
         if ~isempty(strcmp('nli',toPenalise)) && ~isempty(self.c)
            c                 = self.c;
            penalisationStr   = [penalisationStr, '+sigma * max(max(c(x)),0)'];
            n                 = n+1;
         end
         %
         if ~isempty(strcmp('nle',toPenalise)) && ~isempty(self.ceq)
            ceq               = self.ceq;
            penalisationStr   = [penalisationStr, '+sigma * norm(ceq(x),2)^2'];
            n                 = n+1;
         end
         %
         if ~isempty(strcmp('bounds',toPenalise)) && any(self.xMin > -inf) || any(self.xMax < inf)
            xMin              = self.xMin;
            xMax              = self.xMax;
            penalisationStr   = [penalisationStr,'+ sigma * max(max(xMin - x),0)'];
            penalisationStr   = [penalisationStr,'+ sigma * max(max(x - xMax),0)'];
            n                 = n+1;
         end
         fun   = eval(['@(x,sigma) 0',penalisationStr]);
      end
      %% list
      function out = list(self)
         % list the constraints that are present in the constraint object
         out = {};
         if ~isempty(self.A)
            out{end+1} = 'li';
         end
         if ~isempty(self.Aeq)
            out{end+1} = 'le';
         end
         if any(self.xMin > -inf) || any(self.xMax < inf)
            out{end+1} = 'bounds';
         end
         if ~isempty(self.c)
            out{end+1} = 'nli';
         end
         if ~isempty(self.ceq)
            out{end+1} = 'nle';
         end
      end
      %% projectOnBounds
      function X = projectOnBounds(self, X, m ,M)
         % project the population X on the bounds m and M (or by default
         % the bounds of the constraint object)
         if nargin < 3
            m = self.xMin;
            M = self.xMax;
         end
         for i = 1:size(X,2)
            X(:,i) = min(M, max(m,X(:,i)));
         end
      end
      %% adjustBoundsDimension
      function adjustBoundsDimension(self, xDim)
         % adjust the dimension of the bounds to that it matches the actual
         % dimension of x (not known when the constraint object is created)
         if xDim == 1
            return
         end
         if isscalar(self.xMin)
            self.xMin = self.xMin * ones(xDim, 1);
         end
         if isscalar(self.xMax)
            self.xMax = self.xMax * ones(xDim, 1);
         end
      end
      %% getActiveLIC
      function out = getActiveLIC(self, x)
         % returns the active linear inequality constraints
         out = self.A * x > self.b;
      end
      %% getLICFeasibleDir
      function [D,idx] = getLICFeasibleDir(self, X)
         % returns feasible direction for each point of the population X.
         % Returns also a vector of boolean indicating by false the points
         % for which any direction is feasible.
         D     = zeros(size(X));
         idx   = zeros(size(X,2),1);
         for i = 1:size(X,2)
            active = self.getActiveLIC(X(:,i));
            if any(active)
               idx(i) = 1;
               D(:,i) = -1/sum(active) * sum(self.A(active,:))';               
            end
         end
      end
      %% licViolation
      function mv = licViolation(self, x, enableSkip)
         % returns the max violation of the linear inequality constraints
         mv = 0;
         if nargin < 3
            enableSkip = true;
         end
         if enableSkip && any(strcmp('li',self.skip))
            return
         end
         if ~isempty(self.A)
            mv = max(max(self.A * x - self.b),0);
         end
      end
      %% lecViolation
      function mv = lecViolation(self, x, enableSkip)
         % returns the max violation of the linear equality constraints
         mv = 0;
         if nargin < 3
            enableSkip = true;
         end
         if enableSkip && any(strcmp('le',self.skip))
            return
         end
         if ~isempty(self.Aeq)
            mv = norm(self.Aeq * x - self.beq,inf);
         end
      end
      %% nlicViolation
      function mv = nlicViolation(self, x, enableSkip)
         % returns the max violation of the non-linear inequality
         % constraints
         mv = 0;
         if nargin < 3
            enableSkip = true;
         end
         if enableSkip && any(strcmp('nli',self.skip))
            return
         end
         if ~isempty(self.c)
            mv = max(max(self.c(x)),0);
         end
      end
      %% nlecViolation
      function mv = nlecViolation(self, x, enableSkip)
         % returns the max violation of the non-linear equality constraints
         mv = 0;
         if nargin < 3
            enableSkip = true;
         end
         if enableSkip && any(strcmp('nle',self.skip))
            return
         end
         if ~isempty(self.ceq)
            mv = norm(self.ceq(x),inf);
         end
      end
      %% boundViolation
      function mv = boundViolation(self, x, enableSkip)
         % returns the max violation of the bound constraints
         mv = 0;
         if nargin < 3
            enableSkip = true;
         end
         if enableSkip && any(strcmp('bounds',self.skip))
            return
         end
         % min
         xm = self.xMin;
         if isscalar(xm)
            xm = xm * ones(size(x));
         end
         mv = max(mv, max(max(-x + xm),0));
         % max
         xM = self.xMax;
         if isscalar(xM)
            xM = xM * ones(size(x));
         end
         mv = max(mv, max(max(x-xM),0));
      end
      %% containsBox
      function out = containsBox(self)
         % test whether the bounds form a closed box
         out = all(self.xMax < inf) && all(self.xMin > -inf);
      end
      %% maxViolation
      function mV = maxViolation(self, X, enableSkip)
         % returns the max violation of the constraints
         if nargin < 3
            enableSkip = true;
         end
         mV = zeros(1,size(X,2));
         for i = 1:size(X,2)
            v     = -inf;
            x     = X(:,i);
            % Linear inequality constraints
            v     = max(v, self.licViolation(x, enableSkip));
            % Linear equality constraints
            v     = max(v, self.lecViolation(x, enableSkip));           
            % Nonlinear inequality constraints
            v     = max(v, self.nlicViolation(x, enableSkip));
            % Nonlinear equality constraints
            v     = max(v, self.nlecViolation(x, enableSkip));
            % Bound constraints
            v     = max(v, self.boundViolation(x, enableSkip));
            %
            mV(i) = v;
         end
         self.lastViolation = mV;
      end
      %% satisfied
      function [out, msg, mv] = satisfied(self, x, enableSkip)
         % test whether the constraints are satisfied up to some tolerance
         if nargin < 3
            enableSkip = true;
         end
         mv    = self.maxViolation(x, enableSkip);
         out   = mv <= self.tolCst;
         msg   = '';
         if ~out
            msg = sprintf('Maximum constraint violation is above tolerance: %1.2e (tol: %1.2e)',mv, self.tolCst);
         end
      end
   end
   
end

