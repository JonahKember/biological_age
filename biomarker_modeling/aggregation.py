import numpy as np
import pandas as pd

from lifelines import CoxPHFitter
from sklearn.impute import KNNImputer
from scipy.interpolate import RegularGridInterpolator


def _get_delta_age_estimator(delta_age_array):

    # Load delta-age array.
    delta_age    = pd.read_csv(delta_age_array, index_col=0)
    deltas       = delta_age.to_numpy()
    feature_vals = delta_age.index.values
    age_vals     = delta_age.columns.values.astype(float)

    # Fit interpolator.
    delta_age_intepolator = RegularGridInterpolator(points=[feature_vals, age_vals], values=deltas)

    # Define function for performing interpolation.
    def delta_age_estimator(feature_val, age_val):

        # Clip age and feature at min.max of defined range.
        feature_val = np.clip(feature_val, np.min(feature_vals), np.max(feature_vals))
        age_val = np.clip(age_val, np.min(age_vals), np.max(age_vals))

        # Interpolate delta-age at age-feature pair.
        delta_age = delta_age_intepolator((feature_val, age_val))

        return delta_age

    return delta_age_estimator


def get_delta_age_dataframe(config):

    # Merge all biomarker delta-ages.
    output_dir = config['output_dir']
    features = config['features']

    # Load dataframe including all features and age.
    df = pd.read_csv(config['dataframe'])

    # Impute missing data.
    imputer = KNNImputer()
    df_imputed = pd.DataFrame(imputer.fit_transform(df[features]), columns=features)
    df_imputed['age'] = df['age']
    df_imputed.to_csv(f'{output_dir}/data_imputed.csv', index=False)

    # Calculate delta-age for each feature.
    df_delta_ages = pd.DataFrame()
    for feature in features:
        delta_age_estimator = _get_delta_age_estimator(f'{output_dir}/{feature}/model/delta_age-array.csv')
        df_delta_ages[feature] = delta_age_estimator(df_imputed[feature], df_imputed['age'])

    df_delta_ages.to_csv(f'{output_dir}/biomarker_delta_ages.csv', index=False)

    return df_delta_ages


def get_biological_age(config):

    output_dir = config['output_dir']
    df = pd.read_csv(config['dataframe'])
    df_deltas = pd.read_csv(f'{output_dir}/biomarker_delta_ages.csv')

    # Calculate final biological-age.
    df_cox = pd.DataFrame()
    df_cox['duration'] = df[config['duration_col']]
    df_cox['event'] = df[config['event_col']]
    df_cox['deltas'] = df_deltas.sum(axis=1)

    cph = CoxPHFitter().fit(df_cox, 'duration', 'event')
    hazard_ratio = cph.summary['exp(coef)'].item()
    hazard_ratio_conversion = (np.log(hazard_ratio) / np.log(2)) * config['mrdt']

    delta_age = df_cox['deltas'] * hazard_ratio_conversion
    bio_age = df['age'] + delta_age

    df_deltas['delta_age'] = delta_age
    df_deltas['bio_age'] = bio_age
    df_deltas.to_csv(f'{output_dir}/biomarker_delta_ages.csv', index=False)