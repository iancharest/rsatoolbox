% display_singleton_struct(struct_in)
%
% For a struct struct_in, which is not a struct array, and whose contents
% are either non-array structs or not structs, this will display a nested
% representation of the struct hierarchy.
%
% CW 2015-06
function display_singleton_struct(struct_in, prefix)
    import rsa.*
    import rsa.util.*
    
    % Don't need to specify prefix on first call.
    if ~exist('prefix', 'var')
        prefix = '';
    end
    
    field_names = fields(struct_in);
    
    for field_i = 1:numel(field_names)
        
        field_name = field_names{field_i};
        field_contents = struct_in.(field_name);
        
        prints('%s%s', prefix, field_name);
        
        % Prefixes are used for indentation.
        % We ensure that any further indentation starts from the current
        % indent.
        next_prefix = [prefix, '  '];
        
        if isstruct(field_contents)
            % Recursive call.
            display_singleton_struct(field_contents, next_prefix);
        else
            % Wow, is this really the only way to do this?
            type_agnostic_string = evalc('disp(field_contents)');
            
            % Trim trailing newline:
            type_agnostic_string = type_agnostic_string(1:end-1);
            
            % Trim leading whitespace:
            type_agnostic_string = strtrim(type_agnostic_string);
            
            prints('%s%s', next_prefix,  type_agnostic_string);
        end
        
    end
        
end
