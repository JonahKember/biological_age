import os
import json
import numpy as np
import pandas as pd

from glob import glob


def write_metadata(feature_params):

    data_path = feature_params['input_data']
    feature = feature_params['feature']

    df = pd.read_csv(data_path)
    meta = {
        "feature": feature_params['feature'],
        "units": feature_params['units'],
        "label": feature_params['label'],
        "duration_col": feature_params['duration_col'],
        "event_col": feature_params['event_col'],
        "ped_by": feature_params['ped_by'],
        "n_grid": feature_params['n_grid'],
        "n_obs": len(df),
        "n_deaths": int(df[feature_params['event_col']].sum()),
        "feature_min": float(df[feature].min()),
        "feature_max": float(df[feature].max()),
        "feature_mean": float(df[feature].mean()),
        "feature_std": float(df[feature].std()),
    }

    with open(os.path.join(feature_params['model_dir'], "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    return

