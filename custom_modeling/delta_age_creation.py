import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def write_delta_age_array(config, output_dir, mrdt=8, age_vals=np.linspace(30,60,100)):

    feature_vals = np.linspace(
        config['feature_min'],
        config['feature_max'], 100
    )

    # Get expected value given age.
    expected_values = config['reference_value'] + (config['yearly_change'] * (age_vals - config['reference_age']))
    interp_func = interp1d(
        age_vals,
        expected_values,
        kind='linear', fill_value='extrapolate'
    )

    # Get grid of delta-age value.
    age_grid, feature_grid = np.meshgrid(age_vals, feature_vals)
    expected_grid = interp_func(age_grid)

    doubling_time = np.log(config['unit_HR'])/np.log(2)
    delta_age_grid = doubling_time * mrdt * (feature_grid - expected_grid)

    # Write to .csv
    df = pd.DataFrame(delta_age_grid, index=feature_vals, columns=age_vals)
    df.to_csv(f"{output_dir}/{config['feature']}/model/delta_age-array.csv")

    # Save plot.
    fig, ax = plt.subplots(figsize=(10,6))

    c_lim = np.max(np.abs(delta_age_grid))
    contour = ax.contourf(age_vals, feature_vals, delta_age_grid, levels=100, cmap='coolwarm', vmin=-c_lim, vmax=c_lim)
    ax.contour(age_vals, feature_vals, delta_age_grid, levels=10, colors='w', linewidths=0.2)
    ax.figure.colorbar(contour, ax=ax, label='$\Delta$ age', shrink=0.8)

    ax.plot(age_vals, expected_values, linestyle='--', lw=3, color='k')

    ax.set_xlabel('age')
    ax.set_ylabel(config['feature'])

    fig.savefig(f"{output_dir}/{config['feature']}/plots/delta_age-array.png", dpi=1000, bbox_inches='tight')
    plt.close(fig)
