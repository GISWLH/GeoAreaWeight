function o = gaw_parse_options(defaults, args, fname)
%GAW_PARSE_OPTIONS  Parse name-value pairs (case-insensitive) into a struct. (internal)
    o = defaults;
    names = fieldnames(defaults);
    if mod(numel(args), 2) ~= 0
        error('gaw:option', '%s: options must be name-value pairs', fname);
    end
    for k = 1:2:numel(args)
        key = args{k};
        if ~ischar(key) && ~(isstring(key) && isscalar(key))
            error('gaw:option', '%s: option names must be text', fname);
        end
        hit = find(strcmpi(char(key), names), 1);
        if isempty(hit)
            error('gaw:option', '%s: unknown option ''%s''. Valid options: %s', ...
                  fname, char(key), strjoin(names', ', '));
        end
        o.(names{hit}) = args{k + 1};
    end
end
