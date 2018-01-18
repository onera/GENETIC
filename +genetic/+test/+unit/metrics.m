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
classdef metrics <  matlab.unittest.TestCase
   properties
      tol = 1e-12;
   end
   methods (Test)
      %%
      function l2metric(testCase)
         n = 3;
         X = randn(n,4);
         Y = randn(n,5);
         D = genetic.population.metrics.distances(X, Y, 'l2');
         for i = 1:size(X,2)
            for j = 1:size(Y,2)
               expected = norm(X(:,i) - Y(:,j));
               testCase.verifyEqual(D(i,j),expected,'AbsTol',testCase.tol);
            end
         end
      end
      %%
      function sql2metric(testCase)
         n = 3;
         X = randn(n,4);
         Y = randn(n,5);
         D = genetic.population.metrics.distances(X, Y, 'sql2');
         for i = 1:size(X,2)
            for j = 1:size(Y,2)
               expected = norm(X(:,i) - Y(:,j))^2;
               testCase.verifyEqual(D(i,j),expected,'AbsTol',testCase.tol);
            end
         end
      end
      %%
      function l1metric(testCase)
         n = 3;
         X = randn(n,4);
         Y = randn(n,5);
         D = genetic.population.metrics.distances(X, Y, 'l1');
         for i = 1:size(X,2)
            for j = 1:size(Y,2)
               expected = norm(X(:,i) - Y(:,j),1);
               testCase.verifyEqual(D(i,j),expected,'AbsTol',testCase.tol);
            end
         end
      end
      %%
      function linfmetric(testCase)
         n = 3;
         X = randn(n,4);
         Y = randn(n,5);
         D = genetic.population.metrics.distances(X, Y, 'linf');
         for i = 1:size(X,2)
            for j = 1:size(Y,2)
               expected = norm(X(:,i) - Y(:,j),inf);
               testCase.verifyEqual(D(i,j),expected,'AbsTol',testCase.tol);
            end
         end
      end
      %%
      function cosinemetric(testCase)
         n = 3;
         X = randn(n,4);
         Y = randn(n,5);
         D = genetic.population.metrics.distances(X, Y, 'cos');
         for i = 1:size(X,2)
            for j = 1:size(Y,2)
               expected = 1 - X(:,i)'*Y(:,j)/(norm(X(:,i))*norm(Y(:,j)));
               testCase.verifyEqual(D(i,j),expected,'AbsTol',testCase.tol);
            end
         end
      end
   end
end