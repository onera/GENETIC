% CMAES - Covariance Matrix Adaptation Evolution Strategy
%
%  CMA-ES minimizes any objective function by sampling candidate solutions
%  within an adapted  multivariate normal distribution. In order to adapt
%  iteratively the parameters of the search distribution, this algorithm
%  relies on two main principles:
%
%  1/ a maximum-likelihood principle: the mean of the normal distribution
%  is updated s.t. the likelihood of the previously successful candidate
%  solutions is maximized. The covariance matrix of the distribution is
%  updated incrementally such that the likelihood of previously successful
%  search steps is increased. Both updates can be interpreted as a natural
%  gradient descent. Consequently, CMA-ES conducts an iterated principal
%  components analysis of successful search steps while retaining all the
%  principal axes.
%
%  2/ a two-paths evolution principle: two paths of the time evolution of
%  the normal distribution are recorded. These latter contain significant
%  information about the correlation between consecutive steps. Namely, if
%  consecutive steps are taken in a similar direction, the evolution paths
%  become long. The evolution paths are exploited in two ways. One path is
%  used for the covariance matrix adaptation procedure in place of single
%  successful search steps and facilitates a possibly much faster variance
%  increase of favorable directions. The other path is used to conduct an
%  additional step-size control. This latter path aims to make consecutive
%  movements of the distribution mean orthogonal in expectation. The step
%  size control effectively prevents premature convergence yet allowing fast
%  convergence to an optimum.
%
% References
%  [1] A. Ostermeier, A. Gawelczyk, and N. Hansen, A derandomized approach
%      to self-adaptation of evolution strategies. Evolutionary Computation,
%      1994.
%  [2] N. Hansen, S.D. Muller, and P. Koumoutsakos, Reducing the Time
%      Complexity of the Derandomized Evolution Strategy with Covariance Matrix
%      Adaptation (CMA-ES), Evolutionary Computation, 2003.
%  [3] C. Igel, N. Hansen, and S. Roth, Covariance Matrix Adaptation for
%      Multi-objective Optimization, Evolutionary Computation, 2007.
%

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
classdef cmaes < genetic.optimizer.mono & genetic.optimizer.simpleScheme
   properties
      % Control parameter
      nUpdates          = 0;
      % Selection strategy parameters
      offspringSize
      recombinationSize
      recombinationWgts
      varEffectiveness
      % Adaptation strategy parameters
      covCumulation
      stdCumulation
      learningR1Cov
      learningMuCov
      stdDampFactor
      % Dynamical strategy parameters
      evoPathCov
      evoPathStd
      coordinateSystem
      stdMatrix
      covMatrix
      stdInvMatrix
      % Current and previous mean points
      currentMean
      previousMean
      % Path step size
      stepSize
      expNorm
      ga1Stretching      = 7/2;
      ga2Stretching      = 1/2;
      mu0Stretching      = .5e-2;
      %
   end
   methods
      %% constructor
      function self = cmaes(varargin)
         self@genetic.optimizer.mono(varargin{:});
         self@genetic.optimizer.simpleScheme();
         self.methodName = 'cmaes';
         self.longName     = 'Covariance Matrix Adaptation Evolution Strategy';

         self.needFiniteBounds   = true;
         self.initAllAttributes(self.xDim, self.popSize)
         
      end
      %% evolve
      % Return a new group by applying specific mechanisms (proper to this optimizer)
      function group = evolve(obj, group)
         popSize = group.len;
         % empty flag for display
         % increase iteration number
         % update previous mean point
         obj.previousMean  = obj.currentMean;
         % sampled objective value(s) sorting
         tmp               = group.getObjective();
         [~, indexes]      = sort(tmp);
         % sampled point(s) coordinate(s)
         sampledPts        = group.getValue();
         % new mean point calculation
         obj.currentMean   = sampledPts(:,indexes(1:obj.recombinationSize))*obj.recombinationWgts;
         % normalized mean point displacement
         meanVariation     = (obj.currentMean - obj.previousMean)/obj.stepSize;
         %
         ineStd            = obj.stdCumulation*obj.evoPathStd;
         corStd            = sqrt((1 - obj.stdCumulation^2)*obj.varEffectiveness);
         % standard deviation path evolution
         obj.evoPathStd    = ineStd + corStd*obj.stdInvMatrix*meanVariation;
         % standard deviation path norm
         normPathStd       = norm(obj.evoPathStd);
         %
         ineCov            = obj.covCumulation*obj.evoPathCov;
         corCov            = sqrt((1 - obj.covCumulation^2)*obj.varEffectiveness);
         %pow               = 2*group.nEvalGroup/popSize;
         pow               = 2*obj.nEval/popSize;
         
         den               = obj.expNorm;
         %
         tmp               = (normPathStd/(sqrt(1 - obj.stdCumulation^pow)*den)) < (1.4 + 2/(group.xDim + 1));
         % covariance path evolution
         obj.evoPathCov    = ineCov + tmp*corCov*meanVariation;
         % covariance path norm
         norm2PathCov      = obj.evoPathCov*obj.evoPathCov';
         % normalized vectors between sampled point(s) and previous mean point
         delta             = sampledPts(:,indexes(1:obj.recombinationSize)) - repmat(obj.previousMean,1,obj.recombinationSize);
         normalPts         = delta/obj.stepSize;
         %
         ineCovMx          = (1 - obj.learningR1Cov - obj.learningMuCov)*obj.covMatrix;
         corR1CovMx        = norm2PathCov + (1 - tmp)*(1 - obj.covCumulation^2)*obj.covMatrix;
         corMuCovMx        = normalPts*diag(obj.recombinationWgts)*normalPts';
         % covariance matrix adaptation
         obj.covMatrix     = ineCovMx + obj.learningR1Cov*corR1CovMx + obj.learningMuCov*corMuCovMx;
         %
         a                 = (1 - obj.stdCumulation)/obj.stdDampFactor;
         b                 = normPathStd/obj.expNorm;
         % step size update
         obj.stepSize      = obj.stepSize*exp(a*(b - 1));
         %
         diff = obj.nGen - obj.nUpdates;
         cond = diff > popSize/(obj.learningR1Cov + obj.learningMuCov)/group.xDim/10;
         if cond
            obj.nUpdates                  = obj.nGen;
            obj.covMatrix                 = triu(obj.covMatrix) + triu(obj.covMatrix,1)';
            [obj.coordinateSystem,tmp]    = eig(obj.covMatrix);
            obj.stdMatrix                 = sqrt(diag(tmp));
            obj.stdInvMatrix              = obj.coordinateSystem*diag(obj.stdMatrix.^-1)*obj.coordinateSystem';
         end
         % allocate matrix <tmp>:
         tmp = zeros(group.xDim,popSize);
         % sampling (normal law)
         for j = 1:popSize
            tmp(:,j) = obj.currentMean + obj.stepSize*obj.coordinateSystem*(obj.stdMatrix.*randn(group.xDim,1));           % >> m + sig*NormalLaw(0,C) ~ m + sig*sqrt(C)*NormalLaw(0,I) << %
            tmp(:,j) = min(obj.xMax,max(obj.xMin,tmp(:,j)));
         end
         % new sampled points
         group.moveTo(tmp);
      end
      %% initPopulation
      function X0 = initPopulation(self, x0)
         %
         if isempty(x0)
            m = self.currentMean;
         else
            m = mean(x0, 2);
         end
         X0 = zeros(self.xDim, self.popSize);
         for i = 1:self.popSize - size(x0,2)
            X0(:,i)  = m + self.stepSize * self.coordinateSystem * (self.stdMatrix .* randn(self.xDim,1));
         end
         X0             = self.constraints.projectOnBounds(X0);
         if ~isempty(x0)
            X0(:,i+1:self.popSize)  = x0;
         end
      end
      %%
      function initAllAttributes(obj, xDim, popSize)
         obj.currentMean        = obj.xMin + (obj.xMax - obj.xMin).*rand(xDim,1);
         obj.previousMean       = obj.currentMean;
         obj.stepSize           = mean((obj.xMax - obj.xMin)/2);
         obj.evoPathCov         = zeros(xDim,1);
         obj.evoPathStd         = zeros(xDim,1);
         obj.coordinateSystem   = eye(xDim);
         obj.stdMatrix          = ones(xDim,1);
         obj.covMatrix          = obj.coordinateSystem*diag(obj.stdMatrix.^2)*obj.coordinateSystem';
         obj.stdInvMatrix       = obj.coordinateSystem*diag(obj.stdMatrix.^-1)*obj.coordinateSystem';
         obj.recombinationSize  = floor(popSize/2);
         obj.recombinationWgts  = log(popSize/2+1/2) - log(1:popSize/2)';
         obj.recombinationWgts  = obj.recombinationWgts/sum(obj.recombinationWgts);
         obj.varEffectiveness   = sum(obj.recombinationWgts)^2/sum(obj.recombinationWgts.^2);
         obj.covCumulation      = 1 - (4 + obj.varEffectiveness/xDim)/(4 + 2*obj.varEffectiveness/xDim + xDim);
         obj.stdCumulation      = 1 - (2 + obj.varEffectiveness)/(5 + obj.varEffectiveness + xDim);
         obj.learningR1Cov      = 2/(obj.varEffectiveness + (1.3 + xDim)^2);
         obj.learningMuCov      = (obj.varEffectiveness + 1/obj.varEffectiveness - 2)/(obj.varEffectiveness + (2 + xDim)^2);
         obj.learningMuCov      = min(1-obj.learningR1Cov,2*obj.learningMuCov);
         obj.stdDampFactor      = 2 + 2*max(0,sqrt((obj.varEffectiveness - 1)/(1 + xDim))-1) - obj.stdCumulation;
         obj.expNorm            = sqrt(xDim)*(1 - 1/4/xDim + 1/21/xDim^2);
      end
   end
end