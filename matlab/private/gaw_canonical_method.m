function m = gaw_canonical_method(method)
%GAW_CANONICAL_METHOD  Map a method name or alias to 'band'|'cos'|'ellipsoid'|'none'. (internal)
    switch lower(char(method))
        case {'band', 'sphere', 'sphere_band'}, m = 'band';
        case 'cos',                              m = 'cos';
        case {'ellipsoid', 'wgs84'},             m = 'ellipsoid';
        case {'none', 'equal', 'unweighted'},    m = 'none';
        otherwise
            error('gaw:method', 'unknown method ''%s''; expected band, cos, ellipsoid or none', char(method));
    end
end
