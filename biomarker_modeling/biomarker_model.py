import os
import pickle
from biomarker_modeling import model_fitting, utils, plotting, delta_age


class BiomarkerModel:
    """
    A biomarker-specific piecewise hazard model derived from PAMM output.
    """

    def __init__(self, feature, duration_col, event_col, data_path, model_dir, n_bootstraps=10, ped_by=1, n_grid=100, feature_info=None):
        """
        Parameters
        ----------
        feature : str
            Name of the biomarker.
        duration_col : str
            Column in data corresponding to follow-up duration.
        event_col : str
            Column in data corresponding to event indicator.
        data_path : str
            Path to CSV with cleaned input data.
        model_dir : str
            Top-level folder where the model will be stored.
        ped_by : int
            Time step for piecewise expansion (PED).
        n_grid : int
            Number of points for feature and time grids.
        """

        # Specify and create directories for model.
        self.model_dir = model_dir
        self.plots_dir = os.path.join(self.model_dir,'plots')
        self.data_dir = os.path.join(self.model_dir,'data')

        os.makedirs(self.plots_dir, exist_ok=True)
        os.makedirs(self.data_dir, exist_ok=True)

        self.feature = feature
        self.duration_col = duration_col
        self.event_col = event_col
        self.data_path = data_path
        self.n_bootstraps = n_bootstraps
        self.ped_by = ped_by
        self.n_grid = n_grid
        self.feature_info = feature_info

        # Fit model.
        model_fitting.fit(self)
        delta_age.write_age_expected_values(self)
        delta_age.write_delta_age(self)

        # Save plots.
        plotting.save_plots(self)

        # Save model.
        utils.write_metadata(self)
        with open(os.path.join(self.model_dir, "model.pkl"), "wb") as f:
            pickle.dump(self, f)
