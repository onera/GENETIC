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
classdef simulator <  matlab.unittest.TestCase
   properties
      snormx
      snormxWithGrad
   end
   methods(TestMethodSetup)
      function createSimulators(testCase)
         nobj  = 2;
         xDim  = 10;
         A     = ones(xDim, xDim);
         b     = ones(xDim,1);
         testCase.snormx         = genetic.simulator(@(x) [norm(x)^2,norm(A*x - b)^2], xDim, nobj);
         %
         function [y,g] = snwg(x, A, b)
            y = [norm(x)^2,norm(A*x-b)^2];
            if nargout  > 1
               g = [2*x';2*(A'*(A*x-b))'];
            end
         end
         testCase.snormxWithGrad          = genetic.simulator(@(x) snwg(x,A,b), xDim, nobj);
         testCase.snormxWithGrad.gradObj  = true;
      end
   end
   methods (Test)
      %%
      function maxFunEval(testCase)
         sim               = testCase.snormx;
         sim.maxFunEval    = 3;
         X                 = rand(sim.xDim,4);
         Y                 = sim.eval(X);
         testCase.verifyLessThan(norm(Y(:,1:3)),inf)
         testCase.verifyEqual(Y(:,4),inf*ones(sim.nobj,1))
      end
      %%
      function forwardFiniteDiff(testCase)
         sim            = testCase.snormx;
         sim.fdScheme   = 'f';
         simwg          = testCase.snormxWithGrad;
         X              = rand(sim.xDim,2);
         [~,G1]         = sim.eval(X);
         [~,Gref]       = simwg.eval(X);
         for i = 1:size(G1,3)
            testCase.verifyLessThan(norm(G1(:,:,i) - Gref(:,:,i)),1e-4);
         end
      end
      %%
      function backwardFiniteDiff(testCase)
         sim            = testCase.snormx;
         sim.fdScheme   = 'b';
         simwg          = testCase.snormxWithGrad;
         X              = rand(sim.xDim,2);
         [~,G1]         = sim.eval(X);
         [~,Gref]       = simwg.eval(X);
         for i = 1:size(G1,3)
            testCase.verifyLessThan(norm(G1(:,:,i) - Gref(:,:,i)),1e-4);
         end
      end
      %%
      function centeredFiniteDiff(testCase)
         sim            = testCase.snormx;
         sim.fdScheme   = 'c';
         simwg          = testCase.snormxWithGrad;
         X              = rand(sim.xDim,2);
         [~,G1]         = sim.eval(X);
         [~,Gref]       = simwg.eval(X);
         for i = 1:size(G1,3)
            testCase.verifyLessThan(norm(G1(:,:,i) - Gref(:,:,i)),1e-4);
         end

         %
      end      
   end
end