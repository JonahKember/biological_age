import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from glob import glob


def save_plots(model):
    _plot_cumulative_hazard_array(model)
    _plot_delta_age(model)
    _plot_hazard_ratio(model)


def _get_feature_label(model):

    # Set y-axis label
    if model.feature_info is None:
        feature_label = model.feature
    else:
        label = model.feature_info.get('label', model.feature)
        units = model.feature_info.get('units', '')
        feature_label = f"{label}\n{units}" if units else label

        return feature_label


def _plot_cumulative_hazard_array(model):

    # Load hazard array
    cumulative_hazard_df = pd.read_csv(f'{model.model_dir}/model/cumulative_hazard-array.csv', index_col=0)

    time_vals = cumulative_hazard_df.columns.astype(np.float64)
    feature_vals = cumulative_hazard_df.index.astype(np.float64)
    cumulative_risk = cumulative_hazard_df.to_numpy()

    # Plot
    fig, ax = plt.subplots(figsize=(9, 5))

    contour = ax.contourf(time_vals, feature_vals, cumulative_risk, levels=1000, cmap='jet')
    ax.contour(time_vals, feature_vals, cumulative_risk, levels=15, colors='w', linewidths=0.2)
    fig.suptitle('Cumulative risk', x=0.25, y=.95, fontsize=14)

    ax.figure.colorbar(contour, ax=ax, label='Hazard', shrink=0.8)
    ax.set_xlabel('Years')
    feature_label = _get_feature_label(model)
    ax.set_ylabel(feature_label)

    feature_obs = pd.read_csv(model.data_path)[model.feature]
    ax.set_ylim([np.min(feature_obs), np.max(feature_obs)])

    # Save.
    fig.savefig(f'{model.model_dir}/plots/cumulative_hazard-array.png', dpi=1000, bbox_inches='tight')
    plt.close(fig)



def _plot_delta_age(model):

    delta_age = pd.read_csv(f'{model.model_dir}/model/delta_age-array.csv', index_col=0)
    feature_values = delta_age.index

    age_expected_values = pd.read_csv(f'{model.model_dir}/model/age_expected_feature_values-vector.csv', index_col=0)
    age_values = age_expected_values.index.values
    expected_values = age_expected_values[model.feature]

    c_lim = np.max(np.abs(delta_age))

    fig, ax = plt.subplots(figsize=(10,6))
    contour = ax.contourf(age_values, feature_values, delta_age, levels=10, cmap='coolwarm', vmin=-c_lim, vmax=c_lim)
    ax.contour(age_values,feature_values, delta_age, levels=10, colors='w', linewidths=0.2)
    ax.set_xlabel('age')
    feature_label = _get_feature_label(model)
    ax.set_ylabel(feature_label)

    ax.figure.colorbar(contour, ax=ax, label='$\Delta$ age', shrink=0.8)
    ax.plot(age_values, expected_values, color='k', lw=3, linestyle='--')
    fig.savefig(f'{model.model_dir}/plots/delta_age.png', dpi=1000, bbox_inches='tight')
    plt.close(fig)


def _plot_hazard_ratio(model):

    # # Get hazard ratio.
    # cumulative_hazard_df = pd.read_csv(f'{model.model_dir}/model/cumulative_hazard-array.csv', index_col=0)
    # feature_vals = cumulative_hazard_df.index

    # reference_value = pd.read_csv(f'{model.model_dir}/data/data_wide.csv')[model.feature].mean()
    # reference_idx = np.argmin(np.abs(reference_value - feature_vals)).item()

    # hazard_max_follow_up = cumulative_hazard_df.iloc[:,-1]
    # hazard_ratio = hazard_max_follow_up / hazard_max_follow_up.iloc[reference_idx]


    # # Get bootstrapped standard error of hazard ratio.
    # bootstraps =glob(f'{model.model_dir}/model/bootstrapped_estimates/cumulative_hazard-array*.csv')
    # bs_hazard_ratios = []

    # for bootstrap in bootstraps:

    #     bs_cumulative_hazard_df = pd.read_csv(bootstrap, index_col=0)
    #     bs_hazard_max_follow_up = bs_cumulative_hazard_df.iloc[:,-1]
    #     bs_hazard_ratio = bs_hazard_max_follow_up / bs_hazard_max_follow_up.iloc[reference_idx]
    #     bs_hazard_ratios.append(bs_hazard_ratio.to_numpy())

    # hazard_ratio_err = np.array(bs_hazard_ratios).std(axis=0)

    hazard_ratio_df = pd.read_csv(f'{model.model_dir}/model/hazard_ratio.csv', index_col='feature')
    feature_vals = hazard_ratio_df.index.astype(float)
    hazard_ratio = hazard_ratio_df['HR']
    hazard_ratio_lower = hazard_ratio_df['HR_lower']
    hazard_ratio_upper = hazard_ratio_df['HR_upper']

    # Plot.
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 12], hspace=0.05)

    ax_hist = fig.add_subplot(gs[0])
    ax_main = fig.add_subplot(gs[1], sharex=ax_hist)

    # Plot hazard ratio ± SE.
    ax_main.plot(feature_vals, hazard_ratio, lw=2, color='k')

    ax_main.grid(linestyle='--', alpha=.5)
    ax_main.set_ylabel('Hazard-ratio')

    ax_main.fill_between(
        feature_vals,
        hazard_ratio_lower,
        hazard_ratio_upper,
        facecolor='tab:gray', alpha=.25
    )
    ax_main.axhline(1, color='k', linestyle='--')


    # Add boxplot showing distribution of feature to top of plot.
    feature_obs = pd.read_csv(model.data_path)[model.feature]

    sns.boxplot(feature_obs, color='w', ax=ax_hist, orient='h', showfliers=False)
    ax_hist.xaxis.set_visible(False)
    ax_hist.yaxis.set_visible(False)

    sns.despine(ax=ax_main, top=True)
    sns.despine(ax=ax_hist,left=True, bottom=True)
    feature_label = _get_feature_label(model)
    ax_main.set_xlabel(feature_label)

    # Save
    fig.savefig(f'{model.model_dir}/plots/hazard_ratio.png', dpi=1000, bbox_inches='tight')
    plt.close(fig)
