import numpy as np
import pandas as pd

from pygam import LinearGAM, s


def write_age_expected_values(feature_params):
    """
    Fit age → expected biomarker value and save vector.
    """

    feature = feature_params['feature']
    model_dir = feature_params['model_dir']

    data_wide = pd.read_csv(f'{model_dir}/data/data_wide.csv')
    X = np.asarray(data_wide['age']).reshape(-1, 1)
    y = np.asarray(data_wide[feature])

    # Smooth function: highly penalized to reduce overfitting
    age_expected_value = LinearGAM(s(0), lam=5000).fit(X, y)

    age_values = np.arange(data_wide['age'].min(), data_wide['age'].max(), 1)
    expected_values = age_expected_value.predict(age_values)

    lower, upper = age_expected_value.prediction_intervals(age_values, width=0.68).T
    standard_deviation = (upper - lower) / 2

    expected_values_df = pd.DataFrame(expected_values, index=age_values, columns=[feature])
    expected_values_df.to_csv(f'{model_dir}/model/age_expected_feature_values-vector.csv')

    standard_deviation_df = pd.DataFrame(standard_deviation, index=age_values, columns=[feature])
    standard_deviation_df.to_csv(f'{model_dir}/model/age_expected_feature_SD-vector.csv')

    return age_expected_value


def write_delta_age(feature_params):
    """
    Compute delta age for each feature value relative to age-expected value
    using the mortality rate doubling time (MRDT) and cumulative hazard.
    """

    feature = feature_params['feature']
    mrdt = feature_params['mrdt']
    model_dir = feature_params['model_dir']

    # Load cumulative hazard at max follow-up
    cumulative_hazard_max = pd.read_csv(f'{model_dir}/model/cumulative_hazard-array.csv', index_col=0).iloc[:, -1]

    # Load age-expected feature values
    age_expected_values = pd.read_csv(f'{model_dir}/model/age_expected_feature_values-vector.csv', index_col=0)

    feature_values = cumulative_hazard_max.index.values
    expected_values = age_expected_values[feature].values
    age_values = age_expected_values.index.values

    # Compute delta age using MRDT
    delta_age_matrix = []
    eps = 1e-9  # avoid log(0)

    for expected_val in expected_values:
        expected_hazard = np.interp(expected_val, feature_values, cumulative_hazard_max.values)
        expected_hazard = max(expected_hazard, eps)

        observed_hazards = np.clip(cumulative_hazard_max.values, eps, None)
        delta_age = (np.log(observed_hazards / expected_hazard) / np.log(2)) * mrdt
        delta_age_matrix.append(delta_age)

    delta_age_df = pd.DataFrame(
        np.array(delta_age_matrix).T,
        columns=age_values,
        index=feature_values
    )
    delta_age_df.to_csv(f'{model_dir}/model/delta_age-array.csv')

    return delta_age_df
