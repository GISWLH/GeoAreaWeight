function w = gaw_lon_widths(lon, lonBounds, allowSingle)
%GAW_LON_WIDTHS  Width [deg] of each longitude cell (column); handles the 0/360 or +-180 wrap. (internal)
    lon = double(lon(:)); n = numel(lon);
    if any(~isfinite(lon)), error('gaw:lon', 'lon contains NaN or infinite values'); end
    if ~isempty(lonBounds)
        ed = gaw_edges(lonBounds, n, 'lonBounds');
        w = abs(ed(:, 2) - ed(:, 1));
        w(w > 180) = 360 - w(w > 180);     % a pair such as [359 1] crosses the wrap point
        return
    end
    if n == 1
        if allowSingle, w = 1; return; end
        error('gaw:lon', 'cannot infer the width of a single longitude; pass lonBounds');
    end
    d = mod(diff(lon) + 180, 360) - 180;   % unwrap 360-degree jumps
    lon = [lon(1); lon(1) + cumsum(d)];
    w = abs(diff(gaw_infer_bounds(lon)));
end
