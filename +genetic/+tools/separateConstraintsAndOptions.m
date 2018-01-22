function [constraints, options] = separateConstraintsAndOptions(args)
constraints = [];
options     = [];
if length(args) == 1
   % if there is only one additional input argument, it may be the
   % constraints structure or the options structure.
   % To distinguish the two cases, the fields of the argument are tested.
   constraintsFields = {'A','b','Aeq','beq','c','ceq','xMin','xMax'};
   fields            = fieldnames(args{1});
   if ~isempty(intersect(fields,constraintsFields)) && isempty(setdiff(fields, constraintsFields))
      % If the argument contains at least one of the field specific to
      % constraints and no additional field, then it is the constraints
      % structure
      constraints = args{1};
   else
      % otherwise, it is the options structure
      options     = args{1};
   end
elseif length(args) == 2
   % If there are two arguments, then the options structure is the last
   % one.
   constraints = args{1};
   options     = args{2};
elseif length(args)> 2
   error('Too many input arguments')
end
end