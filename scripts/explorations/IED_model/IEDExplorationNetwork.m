classdef IEDExplorationNetwork < SRNNCellTypePairs
    % Only a controlled connectivity replacement; all dynamics are inherited.
    methods
        function obj = IEDExplorationNetwork(varargin)
            obj@SRNNCellTypePairs(varargin{:});
        end
        function replace_connectivity(obj,W)
            assert(obj.is_built && ~obj.has_run,'Build first; replace before run.');
            assert(isequal(size(W),[obj.n obj.n]) && all(isfinite(nonzeros(W))));
            obj.W=sparse(W);
            obj.cached_params=obj.get_params();
        end
        function replace_setpoints(obj,values)
            assert(obj.is_built && ~obj.has_run,'Replace before run.');
            assert(isequal(size(values),[obj.n 1]) && all(isfinite(values)));
            obj.S_c_vec=values;
            obj.cached_params=obj.get_params();
        end
    end
end
