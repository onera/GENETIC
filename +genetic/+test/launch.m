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
function launch(key)
import matlab.unittest.TestSuite;

if nargin < 1
   key = 'all';
end
switch key
   case 'all'
      suite = TestSuite.fromPackage('genetic.test','IncludingSubpackages',true);
   % add cases here to call just a subset of the tests
   otherwise
      suite = TestSuite.fromClass(eval(['?genetic.test.unit.',key]));
end
results  = run(suite);
rt       = table(results);
disp(rt)
end