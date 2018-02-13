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
classdef simulator < handle
   properties
      f_
      nobj        = 1;
      xDim
      %
      gradObj     = false;
      h           = 1e-6;
      fdScheme    = 'f';
      I
      %
      maxFunEval  = inf;
%    end
%    properties (SetAccess = protected)
      nEval = 0;
   end
   methods
      function self = simulator(f, xDim, nobj)
         self.f_     = f;
         self.xDim   = xDim;
         self.nobj   = nobj;
         %
         self.I      = speye(xDim);
      end
      %% eval
      %  Evaluate the simulator (and its derivative) on the group or
      %  population X.
      function [Y, G] = eval(self, X)
         % 
         requireGrad = nargout > 1;
         %
         if isa(X, 'genetic.population.group')
            X = X.getValue();
         end
         %
         nX = size(X,2);
         Y  = zeros(self.nobj, nX);
         %
         if ~requireGrad
            for i = 1:nX
               Y(:,i) = self.f(X(:,i));
            end
            return
         end
         % Gradient is required
         G = zeros(self.nobj, size(X,1), nX);
         if self.gradObj
            for i = 1:nX
               [Y(:,i), G(:,:,i)] = self.f(X(:,i));
            end
            if self.nobj == 1
               G = G(:);
            end
            return
         end
         % Finite differences are used
         switch self.fdScheme
            case 'f'
               F = @(x) self.fforward(x);
            case 'b'
               F = @(x) self.fbackward(x);
            case 'c'
               F = @(x) self.fcentered(x);
         end
         for i = 1:nX
            [Y(:,i), G(:,:,i)] = F(X(:,i));
         end
         if self.nobj == 1
            G = G(:);
         end
         
      end
      %% fforward
      % Function evaluation + forward differentiation
      function [y,g] = fforward(self, x)
         y = self.f(x);
         g = zeros(self.nobj, self.xDim);
         for i = 1:self.xDim
            xp       = x + self.h * self.I(:,i);
            yp       = self.f(xp);
            g(:,i)   = (yp - y)/ self.h;
         end
      end
      %% fbackward
      % Function evaluation + backward differentiation
      function [y,g] = fbackward(self, x)
         y = self.f(x);
         g = zeros(self.nobj, self.xDim);
         for i = 1:self.xDim
            xm       = x - self.h * self.I(:,i);
            ym       = self.f(xm);
            g(:,i)   = (y - ym)/ self.h;
         end
      end
      %% fcentered
      % Function evaluation + backward differentiation
      function [y,g] = fcentered(self, x)
         y = self.f(x);
         g = zeros(self.nobj, self.xDim);
         for i = 1:self.xDim
            xp       = x + self.h * self.I(:,i);
            xm       = x - self.h * self.I(:,i);
            yp       = self.f(xp);
            ym       = self.f(xm);
            g(:,i)   = (yp - ym)/ (2*self.h);
         end
      end      
      %% maxFunEvalReached
      function [out, msg] = maxFunEvalReached(self)
         out = self.nEval >= self.maxFunEval;
         msg = '';
         if out
            msg = sprintf('Maximum number of function evaluations reached (%d)',self.maxFunEval);
         end
      end
   end
   %
   methods (Access=private)
      function [y,g] = f(self,x)
         if self.nEval >= self.maxFunEval
            y = inf * ones(self.nobj, 1);
            g = [];
            return
         end
         self.nEval = self.nEval + 1;
         %
         if nargout < 2
            y = self.f_(x);
            return
         end
         %
         [y,g] = self.f_(x);
         if self.nobj == 1
            g = g.';
         end
      end
      %
   end
   %%
   methods (Static)
      function nEvals = dispatch(nEval, nCore)
         nEvals = zeros(nCore,1);
         tmp = floor(nEval / nCore);
         if tmp > 0
            nEvals = nEvals + tmp;
         end
         
         r  = rem(nEval, nCore);
         if r > 0
            nEvals(1:r) = nEvals(1:r) + 1;
         end
      end
   end
end