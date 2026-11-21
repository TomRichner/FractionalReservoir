function networks = fresp_pick_networks(available, requested)
% FRESP_PICK_NETWORKS The networks a frequency-response figure draws: all the stage ran, or the ones asked for, in the order asked.
%
%   networks = FRESP_PICK_NETWORKS(D.settings.networks, cfg.networks)
%
% Empty `requested` gives every network of the stage in its order. A name
% the stage did not run is an error (fresp_pick_networks:NoSuchNetwork),
% not a silent drop.
%
% See also: fig_freqresp_network, fig_freqresp_plane, select_conditions

if isempty(requested)
    networks = available;
    return
end
networks = reshape(cellstr(requested), 1, []);
missing = setdiff(networks, available);
if ~isempty(missing)
    error('fresp_pick_networks:NoSuchNetwork', 'The stage ran no network %s (it has %s).', ...
        strjoin(missing, ', '), strjoin(available, ', '));
end
end
