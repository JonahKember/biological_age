import json
import argparse
import pandas as pd
import biomarker_modeling


def main(config_path):

    # Get pipeline parameters from config.json.
    with open(config_path,'r') as f:
        config = json.load(f)

    # Fit model for each biomarker.
    df = pd.read_csv(config['dataframe'])
    biomarker_info = pd.read_csv(config['biomarker_info'], index_col='feature')
    config['features'] = biomarker_info.index[biomarker_info.index.isin(df.columns)]

    # Calculate delta-age for each biomarker.
    biomarker_modeling.aggregation.get_delta_age_dataframe(config)
    biomarker_modeling.aggregation.get_biological_age(config)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run biological age pipeline')
    parser.add_argument('config', type=str, help='Path to config.json')
    args = parser.parse_args()

    main(args.config)
