import os
import json
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter

from biomarker_modeling.biomarker_model import BiomarkerModel
from model_aggregating import aggregation

def main():
    # Get pipeline parameters from config.json.
    with open('config.json','r') as f:
        config = json.load(f)

    age_col = config['age_col']
    event_col = config['event_col']
    duration_col = config['duration_col']
    mrdt = config['mrdt']

    output_dir = config['output_dir']
    os.makedirs(output_dir, exist_ok=True)


    # Load data.
    df = pd.read_csv(config['dataframe'])
    df = df.rename(columns={age_col:'age'})
    df = df[df[duration_col] > 0]

    biomarker_info = pd.read_csv('biomarker_info.csv', index_col='feature')
    features = biomarker_info.index
    features = features[features.isin(df.columns)]

    # Fit model for each biomarker.
    for feature in features:

        print(f'Fitting model for {feature}...')

        # Create directory
        model_dir = f'{output_dir}/{feature}/'
        os.makedirs(model_dir, exist_ok=True)

        df_marker = df[[duration_col, event_col, feature, age_col]]

        # Remove values outside 2.5/97.5th percentiles to avoid fitting to extreme outliers.
        df_marker = df_marker[
            (df_marker[feature] >= df_marker[feature].quantile(0.025)) &
            (df_marker[feature] <= df_marker[feature].quantile(0.975))
        ]

        df_marker.to_csv(f'{model_dir}/input_data.csv', index=False)

        model = BiomarkerModel(
            feature=feature,
            duration_col=duration_col,
            event_col=event_col,
            data_path=f'{model_dir}/input_data.csv',
            model_dir=model_dir,
            n_bootstraps=0,
            feature_info={
                'units':biomarker_info.loc[feature,'units'],
                'label':biomarker_info.loc[feature,'label'],
            }
        )


    # Merge all biomarker delta-ages.
    delta_age_paths = [f'{output_dir}/{feature}/model/delta_age-array.csv' for feature in features]
    delta_age_paths = dict(zip(features, delta_age_paths))

    with open(f'{output_dir}/delta_age_paths.json','w') as f:
        json.dump(delta_age_paths, f, indent=4)

    delta_age_df = aggregation.get_delta_age_dataframe(
        config['dataframe'],
        f'{output_dir}/delta_age_paths.json',
        output_dir
    )


    # Calculate final biological-age.
    df_cox = pd.DataFrame()
    df_cox[duration_col] = df[duration_col]
    df_cox[event_col] = df[event_col]
    df_cox['deltas'] = delta_age_df.sum(axis=1)


    cph = CoxPHFitter().fit(df_cox, duration_col, event_col)

    hazard_ratio = cph.summary['exp(coef)'].item()
    hazard_ratio_conversion = (np.log(hazard_ratio) / np.log(2)) * mrdt

    delta_age = df_cox['deltas'] * hazard_ratio_conversion
    bio_age = df['age'] + delta_age

    delta_age_df['delta_age'] = delta_age
    delta_age_df['bio_age'] = bio_age
    delta_age_df.to_csv(f'{output_dir}/biomarker_delta_ages.csv', index=False)


if __name__ == '__main__':
    main()