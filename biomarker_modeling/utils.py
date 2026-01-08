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


# def write_hazard_ratio_SE(model):


#     hazard_ratio_df = pd.read_csv(f'{model.model_dir}/model/cumulative_risk_max_duration-vector.csv', index_col=0)
#     feature_vals = hazard_ratio_df.index

#     # Collect boostrapped estimates of hazard_ratio.
#     bootstrap_files = glob(f'{model.model_dir}/model/bootstrapped_estimates/cumulative_hazard_ratio/*')
#     bootstrap_hazard_ratios = [pd.read_csv(file)['hazard_ratio'].values for file in bootstrap_files]

#     # Calculate boostrapped standard-error.
#     hazard_ratio_SE = np.std(bootstrap_hazard_ratios, axis=0)

#     # Save.
#     hazard_ratio_SE = pd.DataFrame(hazard_ratio_SE, columns=['hazard_ratio_SE'], index=feature_vals)
#     hazard_ratio_SE.to_csv(f'{model.model_dir}/model/bootstrapped_estimates/cumulative_hazard_ratio_SE-vector.csv')

#     return


# def write_optimal_value(model):

#     if model.feature_info is None:
#         units = 'model.feature'
#     else:
#         units = model.feature_info.get('units', '')


#     # Get optimal value (that which minimizes long-term cumulative risk).
#     hazard_ratio_df = pd.read_csv(f'{model.model_dir}/model/cumulative_risk_max_duration-vector.csv', index_col=0)
#     opt_val = hazard_ratio_df['hazard_ratio'].idxmin()

#     # Get a confidence interval for the optimal value using the bootstrapped samples.
#     bootstrap_files = glob(f'{model.model_dir}/model/bootstrapped_estimates/cumulative_hazard_ratio/*')
#     bootstrap_opt_val = [pd.read_csv(file, index_col=0)['hazard_ratio'].idxmin() for file in bootstrap_files]
#     opt_val_SE = np.std(bootstrap_opt_val)

#     # Save the optimal value as pandas series and in human-readable format.
#     pd.Series(
#         [opt_val, opt_val_SE],
#         index=['optimal_value','optimal_value_SE'],
#         name=model.feature
#     ).to_csv(f'{model.model_dir}/model/optimal_value.csv')

#     opt_val_human_readable = f'Optimal-value: {opt_val:.2f} ± {opt_val_SE:.2f} {units}, [{(opt_val-opt_val_SE):.2f}, {(opt_val+opt_val_SE):.2f}]'
#     with open(f'{model.model_dir}/model/optimal_value.txt', 'w') as f:
#         f.write(opt_val_human_readable)

