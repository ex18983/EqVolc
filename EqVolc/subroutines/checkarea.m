function [inarea]=checkarea(lonlims,latlims,eventlon,eventlat)

% Checks whether an event (earthquake, volcano) occurs within the specified
% area limits

% Determine if area boundaries are in [min max] or polygon form
sizebounds=size(lonlims,2);

% Code for [min max] form
if sizebounds==2
    
    % Determine whether area 'wraps around' the boundary at 180E/-180W
    if lonlims(2)>lonlims(1)
        wrap='N';
    elseif lonlims(2)<lonlims(1)
        wrap='Y';
    end
    
    % Determine if in area if area doesn't wrap around
    if wrap=='N'
        if eventlon>=lonlims(1) && eventlon<=lonlims(2) && eventlat>=latlims(1) && eventlat<=latlims(2)
            inarea='Y';
        else
            inarea='N';
        end
        % Determine if in area if area does wrap around (by method of exclusion- ie,
        % seeing if it doesn't lie in the area not covered by the limits)
    elseif wrap=='Y'
        if eventlon>lonlims(2) && eventlon<lonlims (1)
            inarea='N';
        elseif eventlat>=latlims(1) && eventlat<=latlims(2)
            inarea='Y';
        else
            inarea='N';
        end
    end
    
% Code for [polygons] form
elseif sizebounds>2
    
    % Determine if point in question lies within polygon(s)
    in=inpolygon(eventlon,eventlat,lonlims,latlims);
    if in==1
        inarea='Y';
    else
        inarea='N';
    end
    
end

end