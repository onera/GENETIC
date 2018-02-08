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
classdef linesearch < genetic.optimizer.mono
   % Parameters
   %     dir : type of direction to be used. Can be 'steepest', 'bfgs',
   %           'newton' ('bfgs')
   properties
      tolGrad  = 1e-9;
      dir      = 'bfgs';
      % linesearch parameters
      tolAlpha = 1e-14;
      lsMethod = 'strongWolfe';
      % Strong Wolfe LS
      c1       = 0.01;
      c2       = 0.9;
      t1       = 9;
      t2       = 0.1;
      t3       = 0.5
      % backtracking LS
      cBTLS    = 0.5;
      %
   end
   %
   properties(SetAccess = protected)
      %
      alpha0 = 1;
      x0
      y0
      g0
      df0
      alpha_
      LSStop = false;
      %
      d
      %
      nGen = 0;
      iB
      I
   end
   %%
   methods
      function self = linesearch(varargin)
         self@genetic.optimizer.mono(varargin{:});
         %
         self.methodName   = 'linesearch';
         self.longName     = 'Local line search descent method';
         %
         % Adding stopping criterion concerning the linesearch failing to
         % find any point satisfying the conditions
         alphaTolReached = struct('name', 'alphaTolReached', 'cleanStop', true, 'condByConstraints', false);
         self.addStopTest(alphaTolReached);
         % toDisplay
         self.toDisplay{2}.name  = 'iter';
         self.removeDisplay('fGen*');
         %
         self.I   = speye(self.xDim,self.xDim);
         self.iB  = self.I;
         self.popSize = 1;
      end
      %% alphaTolReached
      function [out, msg] = alphaTolReached(self, group)
         out = self.LSStop;
         msg = '';
         if out 
            msg = sprintf('Step size alpha in line search below tolerance: %1.2e (tol: %1.2e)', self.alpha_, self.tolAlpha);
         end
      end
      %%
      function [xopt, fopt, info] = apply(self, x0)
         x0Provided = true;
         if nargin < 2 || isempty(x0)
            x0          = [];
            x0Provided  = false;
         end
         self.printHeader(x0Provided);
         self.printStart();
         % Initialisation -------------------------------------------------
         % Group creation
         group = self.makeGroup();
         % Initial population
         X0       = self.initPopulation(x0);
         group.moveTo(X0);
         % Initial evaluation
         [Y0,G0]  = self.simulator.eval(group);
         group.tellObjective(Y0, G0)
         %
         self.printIterHead()
         self.printIter(group);
         % Main optimization loop -----------------------------------------
         while ~self.stop(group)
            self.nGen      = self.nGen + 1;
            % Descent direction
            x              = group.at(1);
            self.getDirection(x);
            % Line search in the direction d
            self.x0        = x.value;
            self.y0        = x.objective;
            self.g0        = x.grad;
            self.df0       = self.g0' * self.d;
            %
            switch self.lsMethod
               case 'strongWolfe'
                  [self.alpha_, falpha, galpha, errLS]   = self.strongWolfeLS();
                  reevaluate                             = false;
               case 'btls'
                  [self.alpha_,errLS]  = self.backtrackingLS();
                  reevaluate           = true;
            end
            %
            if ~errLS
               xalpha      = self.x0 + self.alpha_ * self.d;
            else
               self.LSStop = true;
               reevaluate  = false;
               xalpha      = self.x0;
               falpha      = self.y0;
               galpha      = self.g0;
            end
            % Moving to the new point
            % xk+1 = xk + alpha * d
            group.moveTo(xalpha)
            %
            if reevaluate
               [falpha, galpha] = self.simulator.eval(group);
            end
            group.tellObjective(falpha, galpha);
            %
            if mod(self.nGen,50) == 0
               self.printIterHead();
            end
            self.printIter(group);
         end
         self.printEnd();
         % Output variable assignement ------------------------------------
         [xopt, fopt] = group.getAbsoluteBest();
         %
         info = self.fillInfo(group);
      end      
      %%
      function getDirection(self, x)
         switch self.dir
            case 'steepest'
               self.d = -x.grad;
            case 'bfgs'
               if ~isempty(x.previousGrad)
                  sk       = x.value - x.previousValue;
                  yk       = x.grad - x.previousGrad;
                  self.iB  = (self.I - sk*yk'/(yk'*sk))*self.iB*(self.I - yk*sk'/(yk'*sk)) + sk*sk'/(yk'*sk);
                  self.d   = -self.iB * x.grad;
               else
                  self.d   = -x.grad;
               end
         end
      end
      %%
      function out = ArmijoTest(self, y)
         out = y > self.y0 + self.c1 * self.alpha_ * self.df0;
      end
      %% strongWolfeLS
      function [alpha, falpha, galpha, err] = strongWolfeLS(self)
         err      = false;
         alpha    = self.alpha0;
         palpha   = 0;
         pfalpha  = self.y0;
         % Phase 1: Admissible interval determination
         while true
            xalpha             = self.x0 + alpha * self.d;
            [falpha, galpha]   = self.simulator.eval(xalpha);
            if falpha > self.y0 + self.c1 * alpha * self.df0 || falpha >= pfalpha
               a           = palpha;
               b           = alpha;
               break;
            else
               dfalpha = self.d'*galpha;
               if abs(dfalpha) <= self.c2 * abs(self.df0)
                  return;
               end
               if dfalpha >=0
                  a           = alpha;
                  b           = palpha;
                  break;
               end
            end
            % new step length
            x1              = 2*alpha - palpha;
            x2              = alpha + self.t1 * (alpha-palpha);
            palpha          = alpha;
            alpha           = self.newAlpha(x1, x2, self.x0, self.d);
            if alpha < self.tolAlpha
               err   = true;
               return
            end
            %
            pfalpha         = falpha;
         end
         % Phase 2: Interval reduction
         % at this point [a,b] is an admissible interval that contains a
         % local minimum in the direction d. This interval must be shrinked
         % to get the new step length alpha
         while 1
            x1              = a + self.t2 * (b-a);
            x2              = b - self.t3 * (b-a);
%             palpha          = alpha;
            alpha           = self.newAlpha(x1, x2, self.x0, self.d);
            if alpha < self.tolAlpha
               err   = true;
               return
            end
            %
            xa             = self.x0 + a*self.d;
            fa             = self.simulator.eval(xa);
            %
            xalpha            = self.x0 + alpha * self.d;
            [falpha, galpha]  = self.simulator.eval(xalpha);
            if falpha > self.y0 + self.c1* alpha * self.df0 || falpha >= fa
               b = alpha;
            else
               dfalpha = self.d'*galpha;
               if abs(dfalpha) <= self.c2 * abs(self.df0)
                  return;
               end
               if (b-a)*dfalpha >=0
                  b = a;
               end
               a = alpha;
            end
            %
         end
      end
      %% newAlpha
      function alpha = newAlpha(self, x1, x2, x0, direction)
         x3 = x2;
         x2 = (x1+x3)/2;
         f1 = self.simulator.eval(x0+x1*direction);
         f2 = self.simulator.eval(x0+x2*direction);
         f3 = self.simulator.eval(x0+x3*direction);
         beta12 = x1^2-x2^2; gamma12 = x1-x2;
         beta23 = x2^2-x3^2; gamma23 = x2-x3;
         beta31 = x3^2-x1^2; gamma31 = x3-x1;
         alpha  = real(1/2*(beta23*f1+ beta31*f2+ beta12*f3)/...
                           (gamma23*f1+gamma31*f2+gamma12*f3));
         
         if isnan(alpha)
            alpha = (x1 + x2)/2;
         end
      end
      %% backtrackingLS
      function [alpha,err] = backtrackingLS(self)
         alpha    = self.alpha0;
         err      = false;
         while true
            xalpha   = self.x0 + alpha * self.d;
            falpha   = self.simulator.eval(xalpha);
            if  falpha < self.y0 + self.c1 * alpha * self.df0;
               break;
            end
            alpha = self.cBTLS * alpha;
            if alpha < self.tolAlpha
               err   = true;
               break;
            end
         end

      end
      
   end
   %%
   methods (Static)
   end
end