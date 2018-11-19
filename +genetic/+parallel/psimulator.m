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
classdef psimulator < genetic.simulator
    properties
        communicator = [];
        nCores = [];
    end
    methods
        function self = psimulator(simIn, com, nCores)
            self.f_                 = simIn.f_;
            self.xDim               = simIn.xDim;
            self.nobj               = simIn.nobj;
            self.I                  = simIn.I;
            self.gradObj            = simIn.gradObj;
            %
            self.communicator       = com;
            self.nCores             = nCores;
        end
        
        function  [Y,G] = eval(self, X)
            %
            G = [];
            requireGrad = nargout > 1;
            %
            if isa(X, 'genetic.population.group')
                X = X.getValue();
            end
            %
            nX = size(X,2);
            Y  = zeros(self.nobj, nX);
            % Gradient is not required
            if ~requireGrad
                idx     = genetic.parallel.dispatch(nX, self.nCores-1);
                data    = cell(length(idx),1);
                for i = 1:length(idx)
                    data{i} = {X(:,idx{i}), idx{i}};
                end
                % Send the data for evaluation to the other nodes
                self.communicator.scatter(1, data); 
                % Get back the evaluated values
                data = self.communicator.gather(1);
                % Assign the values to the output matrix
                for i =1:length(data)
                    Y(:,data{i}{2}) = data{i}{1};
                end
                return
            end
            %% TODO: gradient
        end
        
        function idleBeforeEval(self,communicator)
            if communicator.rank > 1
                % Idle when waiting for evaluation
                while true
                    data    = communicator.receive_src(1); % Receive data from proc 1
                    cmd     = data{1};
                    if cmd == 0
                        ompi.finalize(0);
                    end
                    x   = data{2};
                    idx = data{3};
                    %
                    y   = zeros(self.nobj, length(idx));
                    for j = 1:length(idx)
                        y(:,j)   = self.f_(x(:,j));
                    end
                    communicator.send(1,{y, idx});
                end
            end
        end
    end

end