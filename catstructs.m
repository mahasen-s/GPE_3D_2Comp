function out_struct = catstructs(varargin)
nArgs   = length(varargin);

fields  = {};
vals    = {};

for i=1:nArgs
    fields_i    = fieldnames(varargin{i});
    vals_i      = struct2cell(varargin{i});
    fields      = [fields;fields_i]; 
    vals        = [vals;vals_i];
end
out_struct      = cell2struct(vals,fields);

end