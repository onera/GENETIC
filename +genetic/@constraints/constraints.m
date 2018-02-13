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
      function self = constraints(cstStructure)
         %
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
      %%
      function out = empty(self)
         out = isempty(self.A) && isempty(self.Aeq) && isempty(self.c) && isempty(self.ceq) && all(self.xMin==-inf) && all(self.xMax==inf);
      end
      %% penalise
      function [fun, n] = penalise(self, toPenalise)
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
            penalisationStr   = [penalisationStr,'+ sigma * norm(Aeq * x - beq,inf) '];
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
            penalisationStr   = [penalisationStr, '+sigma * norm(ceq(x),inf)'];
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
         if nargin < 3
            m = self.xMin;
            M = self.xMax;
         end
         for i = 1:size(X,2)
            X(:,i) = min(M, max(m,X(:,i)));
         end
      end
      %%
      function adjustBoundsDimension(self, xDim)
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
      %%
%       function makeBoundsFinite(self, infValue)
%          self.xMin(self.xMin==-inf) = -infValue;
%          self.xMax(self.xMax==inf) = infValue;
%       end
      %%
      function mv = licViolation(self, x, enableSkip)
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
      %%
      function mv = lecViolation(self, x, enableSkip)
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
      %%
      function mv = nlicViolation(self, x, enableSkip)
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
      %%
      function mv = nlecViolation(self, x, enableSkip)
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
      %%
      function mv = boundViolation(self, x, enableSkip)
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
      %%
      function out = containsBox(self)
         out = all(self.xMax < inf) && all(self.xMin > -inf);
      end
      %%
      function mV = maxViolation(self, X, enableSkip)
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
      %%
      function [out, msg, mv] = satisfied(self, x, enableSkip)
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

