import os
import json
import numpy as np
import pandas as pd

from glob import glob


def write_metadata(model):

    df = pd.read_csv(model.data_path)
    meta = {
        "feature": model.feature,
        "duration_col": model.duration_col,
        "event_col": model.event_col,
        "ped_by": model.ped_by,
        "n_grid": model.n_grid,
        "n_obs": len(df),
        "n_deaths": int(df[model.event_col].sum()),
        "feature_min": float(df[model.feature].min()),
        "feature_max": float(df[model.feature].max()),
        "feature_mean": float(df[model.feature].mean()),
        "feature_std": float(df[model.feature].std()),
    }

    if model.feature_info:
        meta["feature_info"] = model.feature_info

    with open(os.path.join(model.model_dir, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    return

