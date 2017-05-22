function struct_out     = appendfields(struct_in,varargin)
fields_append   = varargin(1:2:end)';
vals_append     = varargin(2:2:end)';
struct_append   = cell2struct(vals_append,fields_append);
struct_out      = catstructs(struct_in,struct_append);

end