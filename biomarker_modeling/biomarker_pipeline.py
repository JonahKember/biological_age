import os
from biomarker_modeling import model_fitting, utils, plotting, delta_age


def run_biomarker_model(config, feature):

    # Intialize output directory for current feature.
    model_dir = f"{config['output_dir']}/{feature}"
    os.makedirs(f'{model_dir}/plots', exist_ok=True)
    os.makedirs(f'{model_dir}/data', exist_ok=True)


    model_fitting.fit(config)
    delta_age.write_age_expected_values(config)
    delta_age.write_delta_age(config)

    plotting.save_plots(config)
    utils.write_metadata(config)

