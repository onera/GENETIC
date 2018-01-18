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
classdef constraints <  matlab.unittest.TestCase
   methods (Test)
      %%
      function boundViolation(testCase)
         n     = 3;
         xMin  = zeros(3,1);
         cst   = genetic.constraints(struct('xMin',xMin));
         x     = -1e-3 * ones(n,1);
         mvb   = cst.boundViolation(x);
         mv    = cst.maxViolation(x);
         testCase.verifyEqual(mvb, 1e-3);
         testCase.verifyEqual(mv, 1e-3);
      end
      %%
      function ceqViolation(testCase)
         n     = 3;
         alpha = 1e-3;
         x     = alpha*ones(n,1);
         cst   = genetic.constraints(struct('ceq',@(x) norm(x,inf)));
         mvceq = cst.nlecViolation(x);
         mv    = cst.maxViolation(x);
         testCase.verifyEqual(mvceq, alpha);
         testCase.verifyEqual(mv, alpha);
      end
      %%
      function cViolation(testCase)
         n     = 3;
         alpha = 1e-3;
         x     = alpha*ones(n,1);
         cst   = genetic.constraints(struct('c',@(x) norm(x,inf)));
         mvc   = cst.nlicViolation(x);
         mv    = cst.maxViolation(x);
         testCase.verifyEqual(mvc, alpha);
         testCase.verifyEqual(mv, alpha);
      end
   end
end
