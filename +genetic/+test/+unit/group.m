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
classdef group <  matlab.unittest.TestCase
   %
   methods(Test)
      function extractElements(testCase)
         n = 3;
         N = 4;
         G = genetic.population.group(N);
         X = rand(n, N);
         Y = rand(1,N);
         %
         G.moveTo(X);
         G.tellObjective(Y);
         %
         id          = [1,4,2];
         expectedY   = [];
         for i = 1:length(id)
            expectedY = [expectedY, G.at(id(i)).objective];
         end
         testCase.verifyEqual(expectedY,G.getObjective(id));
         
      end
      %%
      function ringTopology(testCase)
         n = 4;
         T = genetic.population.topology.ring(n);
         testCase.verifyEqual(sort(T{1},'ascend'), [2;4]);
         testCase.verifyEqual(sort(T{2},'ascend'), [1;3]);
         testCase.verifyEqual(sort(T{3},'ascend'), [2;4]);
         testCase.verifyEqual(sort(T{4},'ascend'), [1;3]);
      end
      %%
      function fullTopology(testCase)
         n = 4;
         T = genetic.population.topology.full(n);
         testCase.verifyEqual(sort(T{1},'ascend'), [2;3;4]);
         testCase.verifyEqual(sort(T{2},'ascend'), [1;3;4]);
         testCase.verifyEqual(sort(T{3},'ascend'), [1;2;4]);
         testCase.verifyEqual(sort(T{4},'ascend'), [1;2;3]);
      end
   end
end