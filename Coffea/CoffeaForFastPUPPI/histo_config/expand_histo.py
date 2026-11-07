def expand_histo_yaml(raw: dict) -> dict:
    """
    Expands compact histogram YAML into dict format expected by processor.py.

    Supports:
      defaults:
        <named edges/ranges/etc>
      objects: [pf, puppi, ...]   OR  objects: {pf: {}, puppi: {}}
      plots:
        <plotname>: {variables: [...], axes: [...] , logy?: bool}
      per_object:
        <groupname>:
          objects: [objA, objB, ...]      # NEW: apply these plots to many objects
          plots:
            <plotname>: {variables: [...], axes: [...], logy?: bool}
      special:
        <hname>: {object: <obj>, variables: [...], axes: [...], logy?: bool}

    Also supports axis references:
      edges_ref: <defaults key holding list>
      range_ref: <defaults key holding [min,max]>
    """

    if raw is None:
        raise ValueError("expand_histo_yaml: raw is None (empty YAML?)")

    defaults = raw.get("defaults", {}) or {}
    objects = raw.get("objects", {}) or {}
    plots = raw.get("plots", {}) or {}
    per_object = raw.get("per_object", {}) or {}
    special = raw.get("special", {}) or {}

    # ---- normalize objects ----
    if isinstance(objects, list):
        object_names = list(objects)
    elif isinstance(objects, dict):
        object_names = list(objects.keys())
    else:
        raise TypeError(f"'objects' must be list or dict, got {type(objects)}")

    def _get_logy(cfg):
        return bool(cfg.get("logy", defaults.get("logy", False)))

    def _resolve_axis(ax_cfg):
        ax = dict(ax_cfg)

        if "edges_ref" in ax:
            ref = ax.pop("edges_ref")
            if ref not in defaults:
                raise KeyError(f"edges_ref '{ref}' not found in defaults")
            ax["edges"] = defaults[ref]

        if "range_ref" in ax:
            ref = ax.pop("range_ref")
            if ref not in defaults:
                raise KeyError(f"range_ref '{ref}' not found in defaults")
            ax["range"] = defaults[ref]

        return ax

    def _resolve_axes(axes):
        return [_resolve_axis(a) for a in axes]

    out = {}

    # ---- 1) common plots applied to all objects ----
    for obj_name in object_names:
        for plot_name, cfg in plots.items():
            hname = f"{obj_name}_{plot_name}"
            out[hname] = {
                "object": obj_name,
                "variables": cfg["variables"],
                "axes": _resolve_axes(cfg["axes"]),
                "logy": _get_logy(cfg),
            }

    # ---- 2) per_object extra plots ----
    # allow per_object.<group>.objects: [...]
    for group_name, group_cfg in per_object.items():
        group_cfg = group_cfg or {}
        obj_list = group_cfg.get("objects", None)
        group_plots = group_cfg.get("plots", {}) or {}

        # if no objects list, assume group_name is the object name
        if obj_list is None:
            obj_list = [group_name]

        if not isinstance(obj_list, list):
            raise TypeError(f"per_object.{group_name}.objects must be a list, got {type(obj_list)}")

        for obj_name in obj_list:
            for plot_name, cfg in group_plots.items():
                hname = f"{obj_name}_{plot_name}"
                out[hname] = {
                    "object": obj_name,
                    "variables": cfg["variables"],
                    "axes": _resolve_axes(cfg["axes"]),
                    "logy": _get_logy(cfg),
                }

    # ---- 3) special fully-specified histograms ----
    for hname, cfg in special.items():
        if "object" not in cfg:
            raise KeyError(f"special plot '{hname}' missing required field 'object'")
        if "variables" not in cfg:
            raise KeyError(f"special plot '{hname}' missing required field 'variables'")
        if "axes" not in cfg:
            raise KeyError(f"special plot '{hname}' missing required field 'axes'")

        out[hname] = {
            "object": cfg["object"],
            "variables": cfg["variables"],
            "axes": _resolve_axes(cfg["axes"]),
            "logy": _get_logy(cfg),
        }

    return out
