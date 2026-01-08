import os
import json
import pandas as pd
import biomarker_modeling

def main():

    # Get pipeline parameters from config.json.
    with open('config.json','r') as f:
        config = json.load(f)

    df           = pd.read_csv(config['dataframe'])
    event_col    = config['event_col']
    duration_col = config['duration_col']
    output_dir   = config['output_dir']

    os.makedirs(config['output_dir'], exist_ok=True)
    df = df[df[duration_col] > 0]

    # Fit model for each biomarker.
    biomarker_info = pd.read_csv(config['biomarker_info'], index_col='feature')
    features = biomarker_info.index[biomarker_info.index.isin(df.columns)]

    for feature in features:
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
        df_marker.to_csv(f"{feature_params['model_dir']}/input_data.csv", index=False)
        feature_params['input_data'] = f"{feature_params['model_dir']}/input_data.csv"


        biomarker_modeling.model_fitting.fit(feature_params)
        biomarker_modeling.delta_age.write_age_expected_values(feature_params)
        biomarker_modeling.delta_age.write_delta_age(feature_params)

        biomarker_modeling.plotting.save_plots(feature_params)
        biomarker_modeling.utils.write_metadata(feature_params)





#     # Merge all biomarker delta-ages.
#     delta_age_paths = [f'{output_dir}/{feature}/model/delta_age-array.csv' for feature in features]
#     delta_age_paths = dict(zip(features, delta_age_paths))

#     with open(f'{output_dir}/delta_age_paths.json','w') as f:
#         json.dump(delta_age_paths, f, indent=4)

#     delta_age_df = biomarker_modeling.aggregation.get_delta_age_dataframe(
#         config['dataframe'],
#         f'{output_dir}/delta_age_paths.json',
#         output_dir
#     )


#     # Calculate final biological-age.
#     df_cox = pd.DataFrame()
#     df_cox[duration_col] = df[duration_col]
#     df_cox[event_col] = df[event_col]
#     df_cox['deltas'] = delta_age_df.sum(axis=1)


#     cph = CoxPHFitter().fit(df_cox, duration_col, event_col)

#     hazard_ratio = cph.summary['exp(coef)'].item()
#     hazard_ratio_conversion = (np.log(hazard_ratio) / np.log(2)) * mrdt

#     delta_age = df_cox['deltas'] * hazard_ratio_conversion
#     bio_age = df['age'] + delta_age

#     delta_age_df['delta_age'] = delta_age
#     delta_age_df['bio_age'] = bio_age
#     delta_age_df.to_csv(f'{output_dir}/biomarker_delta_ages.csv', index=False)


if __name__ == '__main__':
    main()