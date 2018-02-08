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
classdef gss < genetic.optimizer.mono & genetic.optimizer.simpleScheme
   % GSS - Generating Set Search Method
   % 
   % The Generating Set Search Method [1,2] is a generic framework for direct 
   % search. It encompasses several well-known schemes such as: coordinate 
   % search (or alternating directional search), Hooke Jeeves pattern 
   % search, etc.
   % Multi-directional search is also a GSS method but it is implemented in
   % mdsearch for convenience.
   %
   % It has been shown that this framework is convergent towards a
   % stationary point of the objective function assuming it is smooth
   % enough.
   %  
   % Parameters
   %     generatingSet  : type of generating set to consider. Can be
   %                      'compass' or 'full' when the dimension is < 7
   %                      ('compass')
   %     concurrentEval : boolean that indicates whether the whole
   %                      generating set should be evaluated at once or
   %                      sequentially (false)
   %     DeltaTol       : tolerance on the length of the step size (1e-8)
   %     Delta          : initial step size (1)
   %     theta          : decrease factor for the step size. It must be < 1
   %                      (1/2)
   %     phi            : increase factor for the step size. It must be > 1
   %                      (2)
   %     rho            : positive function of the step size s.t. a trial
   %                      point xt is deemed acceptable if 
   %                       f(xt) < fbest - rho(Delta)
   %                      It must be a handle function (default is @(x) 0)
   % References
   %     [1] T.G. Kolda, R.M. Lewis and V. Torczon. Optimization by Direct
   %         Search: New Perspectives on Some Classical and Mordern Methods. SIAM
   %         Review. 2003.
   %     [2] V. Torczon. On the convergence of pattern search algorithms. SIAM
   %         J. Optim. 1997.
   %
   properties
      % Method for the generating set
      generatingSet  = 'compass';
      %
      heuristic      = '';
      concurrentEval = false;
      % 
      rho            = @(x) 0;
      phi            = 2;
      theta          = 1/2;
      thetMax        = 0.99;
      Delta          = 1;
      %
      DeltaTol       = 1e-8;
   end
   properties (SetAccess = protected)
      G0             = [];
      KF             = [];
      phik           = [];
      thetak         = [];
      Deltak         = [];
      % Best position
      xkm1           = [];
      fkm1           = [];
      xk             = [];
      fk             = inf;
      % Generating set
      i              = 0;
      nRemDir        = [];
      Gk             = [];
      Hk             = [];
   end
   methods
      %% constructor
      function self = gss(varargin)
         self@genetic.optimizer.mono(varargin{:});
         self@genetic.optimizer.simpleScheme();
         self.methodName   = 'gss';
         self.longName     = 'Generating Set Search';
         self.popSize      = 1;
         %
         switch self.generatingSet
            case 'compass'
               % adsearch
               self.G0  = kron([1,-1],speye(self.xDim,self.xDim));
            case 'full' % For completeness mainly, but should not be used
               % adsearch-like
               if self.xDim >= 6
                  error('The full Generating Set can only be used for very small problems due to its dimension: 3^n')
               end
                  self.G0 = genetic.optimizer.mo.gss.fullGeneratingSet(self.xDim);
            otherwise
               error('unknown')
         end
         self.makeNewDk();
         %
         if self.concurrentEval
            self.popSize = self.lenDk;
         end
         % Adding stopping criterion concerning the size of the simplex
         stepTooSmall = struct('name', 'stepTooSmall', 'cleanStop', true, 'condByConstraints', false);
         self.addStopTest(stepTooSmall);
         % Adding printing function
         stepSize        = struct('name','stepSize','call','printStepSize','align','c','dim',16);
         self.addDisplayElement(stepSize);
         % Setting initial values
         self.phik   = self.phi;
         self.thetak = self.theta;
         self.Deltak = self.Delta;
      end
      function [out, msg] = stepTooSmall(self, group)
         out = self.Deltak < self.DeltaTol;
         msg = '';
         if out
            msg = sprintf('Size of the step below tolerance: %1.2e (tol: %1.2e)',self.Deltak,self.DeltaTol);
         end
         
      end
      %% printStepSize
      function out = printStepSize(self, group)
         out = sprintf('%1.5e',self.Deltak);
      end
      %% Dk
      function d = Dk(self, j)
         d = [];
         if nargin < 2
            d = [self.Hk,self.Gk];
            return
         end
         if j <= size(self.Hk, 2)
            d = self.Hk(:,j);
         elseif j <= size(self.Gk, 2) +size(self.Hk, 2)
            d = self.Gk(:,j - size(self.Hk,2));
         end
      end
      %% lenDk
      function n = lenDk(self)
         n = size(self.Hk,2) + size(self.Gk,2);
      end
      %% nextDirection
      function d = nextDirection(self)
         d = [];
         %
         if ~self.concurrentEval && self.nRemDir > 0
            self.nRemDir   = self.nRemDir -1;
            if self.i+1 > self.lenDk
               idx = rem(self.lenDk, self.i+1);
            else
               idx = self.i+1;
            end
            d              = self.Dk(idx);
            self.i         = self.i+1;
         end
      end
      %% sequentialEval
      function out = sequentialEval(self)
         out = ~self.concurrentEval;
      end
      %% evolve
      function group = evolve(self, group)
         % Trial objective
         fTrial         = group.best.objective;
         success        = fTrial < self.fk - self.rho(self.Deltak);
         %
         dk             = self.nextDirection();
         endOfDkReached = isempty(dk);
         %
         if success
            % Successfull iteration: update current point
            if ~isempty(self.xk)
               self.xkm1 = self.xk;
               self.fkm1 = self.fk;
            end
            self.xk        = group.best.value;
            self.fk        = fTrial;
            % and expand the step size
            self.Deltak    = self.phik * self.Deltak;
            % To restart the set Dk
            endOfDkReached = true;
         else
            % Unsuccessfull iteration
            if endOfDkReached
               % If all the poll of direction has been explored, reduced
               % the step size
               self.Deltak = self.thetak * self.Deltak;
            end
         end
         %
         if endOfDkReached
            Hk = [];
            switch self.heuristic
               case 'HJ'
                  if ~isempty(self.fkm1) && success
                     Hk = 1/self.Deltak * repmat((self.xk - self.xkm1),1,2*self.xDim) + kron([1,-1],speye(self.xDim,self.xDim));
                  end
            end
            %
            self.makeNewDk(Hk)
            %
            dk          = self.nextDirection();
         end
         %
         if self.sequentialEval
            newX  = self.xk + self.Deltak * dk;
         else
            dk    = self.Dk;
            newX  = repmat(self.xk,1,self.lenDk) + self.Deltak * dk;
         end
         group.moveTo(newX);
      end
      %% makeNewDk
      function makeNewDk(self,Hk)
         if nargin < 2
            Hk = [];
         end
         self.Hk        = Hk;
         self.Gk        = self.G0;
         %
         self.nRemDir   = self.lenDk();
         self.i         = 0;
      end

   end
   
   %
   methods (Static)
      %% funnGeneratingSet
      function G = fullGeneratingSet(n)
         Gs    = genetic.optimizer.mo.gss.recursiveGS(n, 1, [], []);
         iter  = 0;
         r     = [];
         c     = [];
         e     = [];
         for i = 1:length(Gs)
            gi    = Gs{i};
            if any(gi{2} ~= 0)
               iter  = iter+1;
               r     = [r;gi{1}];
               c     = [c;iter*ones(length(gi{1}),1)];
               e     = [e;gi{2}];
            end
         end
         G = sparse(r,c,e);
      end
      %% recursiveGS
      function G = recursiveGS(n, j, r, e)
         if j > n
            G = {{r,e}};
            return
         end
         G1    = genetic.optimizer.mo.gss.recursiveGS(n, j+1, [r;j], [e;-1]);
         G2    = genetic.optimizer.mo.gss.recursiveGS(n, j+1, [r;j], [e;0]);
         G3    = genetic.optimizer.mo.gss.recursiveGS(n, j+1, [r;j], [e;1]);
         G     = {G1{:}, G2{:}, G3{:}};
      end
      %
   end
      
end