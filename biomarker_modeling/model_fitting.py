import os
import subprocess


def fit(model):
    """
    Run Rscript fit_model.R and load hazard/cumulative hazard/cumulative risk arrays
    for both time-varying and non-time-varying models.
    """

    script_path = os.path.join(os.path.dirname(__file__), "r", "fit_model.R")

    cmd = [
        "Rscript",
        script_path,
        f"--data_path={model['input_data']}",
        f"--duration_col={model['duration_col']}",
        f"--event_col={model['event_col']}",
        f"--feature={model['feature']}",
        f"--model_dir={model['model_dir']}",
        f"--ped_by={model['ped_by']}",
        f"--n_grid={model['n_grid']}"
    ]

    subprocess.run(cmd, check=True)

