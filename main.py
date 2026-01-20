import os
import json
import argparse
import pandas as pd
import biomarker_modeling

def main(config_path):

    # Get pipeline parameters from config.json.
    with open(config_path,'r') as f:
        config = json.load(f)

    df           = pd.read_csv(config['dataframe'])
    event_col    = config['event_col']
    duration_col = config['duration_col']
    output_dir   = config['output_dir']
    fit_model    = eval(config['fit_model'])

    os.makedirs(config['output_dir'], exist_ok=True)
    df = df[df[duration_col] > 0]

    # Fit model for each biomarker.
    biomarker_info = pd.read_csv(config['biomarker_info'], index_col='feature')
    config['features'] = biomarker_info.index[biomarker_info.index.isin(df.columns)]


    for feature in config['features']:
        if not fit_model: continue

        print(f'Fitting model for {feature}...')

        feature_params = config.copy()
        feature_params['feature'] = feature
        feature_params['units']   = biomarker_info.loc[feature,'units']
        feature_params['label']   = biomarker_info.loc[feature,'label']

        # Intialize output directory for current feature.
        feature_params['model_dir'] = f'{output_dir}/{feature}'
        os.makedirs(f"{feature_params['model_dir']}/plots", exist_ok=True)
        os.makedirs(f"{feature_params['model_dir']}/data", exist_ok=True)

        # Remove values outside 2.5/97.5th percentiles to avoid fitting to extreme outliers.
        df_marker = df[[duration_col, event_col, feature, 'age']]
        df_marker = df_marker[
            (df_marker[feature] >= df_marker[feature].quantile(0.025)) &
            (df_marker[feature] <= df_marker[feature].quantile(0.975))
        ]
        feature_params['input_data'] = f'{output_dir}/{feature}/input_data.csv'
        df_marker.to_csv(feature_params['input_data'], index=False)


        biomarker_modeling.model_fitting.fit(feature_params)
        biomarker_modeling.delta_age.write_age_expected_values(feature_params)
        biomarker_modeling.delta_age.write_delta_age(feature_params)

        biomarker_modeling.plotting.save_plots(feature_params)
        biomarker_modeling.utils.write_metadata(feature_params)

    # Calculate delta-age for each biomarker.
    biomarker_modeling.aggregation.get_delta_age_dataframe(config)
    biomarker_modeling.aggregation.get_biological_age(config)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run biological age pipeline')
    parser.add_argument('config', type=str, help='Path to config.json')
    args = parser.parse_args()

    main(args.config)