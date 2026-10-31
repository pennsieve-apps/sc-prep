"""
Seurat .rds -> counts .h5ad. Run by readers.read_rds in a child process, because the file is
untrusted: `python -I rds.py <in.rds> <out.h5ad> <memory limit in bytes>` with an empty
environment and a timeout. The file is only parsed (rdata, a pure-Python reader of R's
serialization format; nothing in it is evaluated, unlike readRDS in R), and only the slots
needed are converted: the counts of the RNA assay (else the active assay) and meta.data.
Prints one JSON line describing what was read.

Seurat v5 (Assay5): every `counts` layer (`counts`, or `counts.<x>` after split()), placed with
the assay's cells / features maps. Seurat v3/v4 (Assay): the `counts` slot (dgCMatrix).
"""
from __future__ import annotations

import json
import resource
import sys


def main(src: str, dst: str, limit: int) -> None:
    if sys.platform == "linux":         # the image; macOS (local development) does not enforce RLIMIT_AS
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from rdata import parser
    from rdata.conversion import SimpleConverter

    T = parser.RObjectType
    convert = SimpleConverter().convert

    def name_of(tag) -> str:
        t = tag.referenced_object if tag.info.type == T.REF else tag
        return t.value.value.decode()

    def slots(obj) -> dict:
        """An S4 object's slots (its attribute pairlist), unconverted."""
        out, a = {}, obj.attributes
        while a is not None and a.info.type == T.LIST:
            out[name_of(a.tag)] = a.value[0]
            a = a.value[1]
        return out

    def named(obj) -> dict:
        """An R named list, its elements unconverted."""
        names = [str(n) for n in convert(slots(obj)["names"])] if obj.attributes is not None else []
        return dict(zip(names, obj.value))

    def r_class(obj) -> list[str]:
        cls = slots(obj).get("class") if obj.attributes is not None else None
        return [str(c) for c in convert(cls)] if cls is not None else []

    def dgc(obj):
        if "dgCMatrix" not in r_class(obj):
            raise ValueError(f"counts are {'/'.join(r_class(obj)) or obj.info.type.name}, not a dgCMatrix")
        s = slots(obj)
        dim = [int(d) for d in convert(s["Dim"])]
        x = np.asarray(convert(s["x"]))
        return sp.csc_matrix((x, np.asarray(convert(s["i"])), np.asarray(convert(s["p"]))), shape=dim), s

    root = parser.parse_file(src).object
    if root.info.type != T.S4 or "Seurat" not in r_class(root):
        raise ValueError(f"not a Seurat object (class {'/'.join(r_class(root)) or root.info.type.name})")
    top = slots(root)
    assays = named(top["assays"])
    active = str(convert(top["active.assay"])[0])
    name = "RNA" if "RNA" in assays else active
    assay = slots(assays[name])
    assay_class = r_class(assays[name])
    meta = convert(top["meta.data"])

    if "layers" in assay:                                   # Seurat v5 Assay5
        layers = named(assay["layers"])
        keys = [k for k in layers if k == "counts" or k.startswith("counts.")]
        if not keys:
            raise ValueError(f"assay {name} has no counts layer (layers: {', '.join(layers)})")
        cells_map, feats_map = convert(assay["cells"]), convert(assay["features"])
        all_cells = np.asarray(cells_map.coords["dim_0"].values).astype(str)
        all_feats = np.asarray(feats_map.coords["dim_0"].values).astype(str)
        blocks, cells = [], []
        for k in keys:
            m, _ = dgc(layers[k])
            f_idx = np.flatnonzero(np.asarray(feats_map.sel(dim_1=k).values, dtype=bool))
            c_idx = np.flatnonzero(np.asarray(cells_map.sel(dim_1=k).values, dtype=bool))
            if m.shape != (len(f_idx), len(c_idx)):
                raise ValueError(f"layer {k} is {m.shape} but its maps list {len(f_idx)} features x {len(c_idx)} cells")
            t = m.T.tocsr()
            t.indices = f_idx[t.indices].astype(t.indices.dtype)
            blocks.append(sp.csr_matrix((t.data, t.indices, t.indptr), shape=(len(c_idx), len(all_feats))))
            cells.append(all_cells[c_idx])
        X = sp.vstack(blocks, format="csr")
        cells = np.concatenate(cells)
        features = all_feats
    else:                                                   # Seurat v3 / v4 Assay
        keys = ["counts"]
        m, s = dgc(assay["counts"])
        dimnames = convert(s["Dimnames"])
        if m.shape[0] == 0 or dimnames[0] is None or dimnames[1] is None:
            raise ValueError(f"assay {name} has no counts")
        X = m.T.tocsr()
        features, cells = np.asarray(dimnames[0]).astype(str), np.asarray(dimnames[1]).astype(str)
    if len(set(cells)) != len(cells):
        raise ValueError("cell names are not unique across the counts layers")

    meta = meta if isinstance(meta, pd.DataFrame) else pd.DataFrame(index=cells)
    meta.index = meta.index.astype(str)
    meta.columns = [str(c) for c in meta.columns]
    obs = meta.reindex(cells)
    for col in obs.columns:
        if not pd.api.types.is_numeric_dtype(obs[col]) or pd.api.types.is_bool_dtype(obs[col]):
            obs[col] = obs[col].astype("string").fillna("").astype(str)
    import anndata as ad

    ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(features))).write_h5ad(dst)
    print(json.dumps({"assay": name, "assay_class": "/".join(assay_class), "active_assay": active,
                      "layers": keys, "n_cells": int(X.shape[0]), "n_features": int(X.shape[1])}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
