"""
plotter.py
==========
Functions to make DM plots
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

__all__ = [
    'plot_thermal_relic',
    'plot_sigmav_ULs_from_hdf5',
    'plot_sigmav_CLbands_from_hdf5',
    'plot_sigmav_CLbands_as_lines_from_hdf5'
]


def plot_thermal_relic(ax=None, **kwargs):
    """
    Plot the thermal relic, which was taken from Steigman G., Dasgupta B, and Beacom J. F.,
    Precise relic WIMP abundance and its impact onsearches for dark matter annihilation,
    Phys.Rev. D86(2012) 023506, [arXiv:1204.3622].

    Parameters
    ----------
    ax: `matplotlib.pyplot.axes`
    kwargs: kwargs for `matplotlib.pyplot.plot`

    Returns
    -------
    ax: `matplotlib.pyplot.axes`
    """

    ax = plt.gca() if ax is None else ax

    ax.set_xscale('log')
    ax.set_xlabel(r'$m_{\chi} \: [GeV]$')

    ax.set_yscale('log')
    ax.set_ylabel(r'$\langle\sigma v\rangle \, [cm^{3}/s]$')
    
    if 'label' not in kwargs:
        kwargs['label'] =  r'Thermal relic $\langle\sigma v\rangle$'
    
    if 'linestyle' not in kwargs:
        kwargs['linestyle'] = '--'

    if 'linewidth' not in kwargs:
        kwargs['linewidth'] = 1.5
    
    if 'color' not in kwargs:
        kwargs['color'] = 'r'

    # Plot the thermal relic, which was taken from Steigman G., Dasgupta B, and Beacom J. F.,
    # Precise relic WIMP abundance and its impact onsearches for dark matter annihilation,
    # Phys.Rev. D86(2012) 023506, [arXiv:1204.3622]
    thermal_relic_mass = [1.00e-01, 1.78e-01, 3.16e-01, 5.62e-01, 1.00e+00, 1.78e+00, 3.16e+00, 5.62e+00, 1.00e+01, 1.78e+01, 3.16e+01, 5.62e+01, 1.00e+02, 1.78e+02, 3.16e+02, 5.62e+02, 1.00e+03, 1.78e+03,3.16e+03, 5.62e+03, 1.00e+04, 1.00e+05]
    thermal_relic_sigmav = [4.8e-26, 4.9e-26, 5.1e-26, 5.0e-26, 4.7e-26, 4.5e-26, 3.9e-26, 2.8e-26, 2.5e-26, 2.3e-26, 2.2e-26, 2.2e-26, 2.2e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26, 2.3e-26,2.4e-26, 2.4e-26]
    ax.plot(thermal_relic_mass, thermal_relic_sigmav, **kwargs)

    return ax


def plot_sigmav_ULs_from_hdf5(channel,
                            file,
                            key='sigmavULs_Jnuisance',
                            ax=None,
                            **kwargs):
    """
    Plot the sigmav upper limits from hdf5 file.

    Parameters
    ----------
    channel: str
        name of the channel.
    file: `path`
        path to a panda readable hdf5 file.
    key: str
        name of the table of the hdf5 file. 
    ax: `matplotlib.pyplot.axes`
    kwargs: kwargs for `matplotlib.pyplot.plot`

    Returns
    -------
    ax: `matplotlib.pyplot.axes`
    """
    
    ax = plt.gca() if ax is None else ax

    ax.set_xscale('log')
    ax.set_xlabel(r'$m_{\chi} \: [GeV]$')

    ax.set_yscale('log')
    ax.set_ylabel(r'$\langle\sigma v\rangle \, [cm^{3}/s]$')
    
    sigmavULs = pd.read_hdf(file, key='{}/{}'.format(channel, key), mode='r')
    masses = np.squeeze(sigmavULs[['masses']].to_numpy())
    data = np.squeeze(sigmavULs[['data'.format(channel)]].to_numpy())
    
    ax.plot(masses, data, **kwargs)
    
    return ax

def plot_sigmav_CLbands_from_hdf5(channel,
                                file,
                                key='sigmavULs_Jnuisance',
                                ax=None):
    """
    Plot the sigmav confidence limit bands from hdf5 file.

    Parameters
    ----------
    channel: str
        name of the channel.
    file: `path`                   
        path to a panda readable hdf5 file.
    key: str
        name of the table of the hdf5 file.
    ax: `matplotlib.pyplot.axes`

    Returns
    -------
    ax: `matplotlib.pyplot.axes`
    """
    
    ax = plt.gca() if ax is None else ax

    ax.set_xscale('log')
    ax.set_xlabel(r'$m_{\chi} \: [GeV]$')

    ax.set_yscale('log')
    ax.set_ylabel(r'$\langle\sigma v\rangle \, [cm^{3}/s]$')

    sigmavULs = pd.read_hdf(file, key='{}/{}'.format(channel, key), mode='r')
    masses = np.squeeze(sigmavULs[['masses']].to_numpy())
    simulations = len(sigmavULs.columns)-2
    # Calculate the median (null hypothesis)
    simu_tab = sigmavULs.drop(['data','masses'], axis=1).to_numpy()
    null_hypothesis = np.median(simu_tab, axis=1)
    sv_plus1 = np.percentile(simu_tab, 84.14, axis=1)
    sv_plus2 = np.percentile(simu_tab, 97.725, axis=1)
    sv_minus1 = np.percentile(simu_tab, 15.87, axis=1)
    sv_minus2 = np.percentile(simu_tab, 2.275, axis=1)

    ax.plot(masses, null_hypothesis, label=r'$ H_{0} $ median', c='k', linewidth=0.75, linestyle='--')
    ax.fill_between(masses,sv_plus1,sv_minus1, color='green', alpha=0.5, linewidth=0)
    ax.fill_between(masses,sv_plus1,sv_plus2, color='yellow', alpha=0.5, linewidth=0)
    ax.fill_between(masses,sv_minus1,sv_minus2, color='yellow', alpha=0.5, linewidth=0)

    # Creating dummy data, which is needed for beautiful legend
    dummy_val = np.ones(len(masses))*1e-32
    ax.plot(masses,dummy_val,label='$ H_{0} \, 68\% $ containment', c='green', alpha=0.5, linewidth=6)
    ax.plot(masses,dummy_val,label=r'$ H_{0} \, 95\% $ containment', c='yellow', alpha=0.5, linewidth=6)

    return ax



def plot_sigmav_CLbands_as_lines_from_hdf5(channel,
                                                file,
                                                key='sigmavULs_Jnuisance',
                                                ax=None,
                                                **kwargs):
    """
    Plot the sigmav confidence limit bands as lines from hdf5 file.

    Parameters
    ----------
    channel: str
        name of the channel.
    file: `path`
        path to a panda readable hdf5 file.
    key: str
        name of the table of the hdf5 file.
    ax: `matplotlib.pyplot.axes`
    kwargs: kwargs for `matplotlib.pyplot.plot`

    Returns
    -------
    ax: `matplotlib.pyplot.axes`
    """

    ax = plt.gca() if ax is None else ax

    ax.set_xscale('log')
    ax.set_xlabel(r'$m_{\chi} \: [GeV]$')

    ax.set_yscale('log')
    ax.set_ylabel(r'$\langle\sigma v\rangle \, [cm^{3}/s]$')

    sigmavULs = pd.read_hdf(file, key='{}/{}'.format(channel, key), mode='r')
    masses = np.squeeze(sigmavULs[['masses']].to_numpy())
    simulations = len(sigmavULs.columns)-2
    # Calculate the median (null hypothesis)
    simu_tab = sigmavULs.drop(['data','masses'], axis=1).to_numpy()
    null_hypothesis = np.median(simu_tab, axis=1)
    sv_plus1 = np.percentile(simu_tab, 84.14, axis=1)
    sv_plus2 = np.percentile(simu_tab, 97.725, axis=1)
    sv_minus1 = np.percentile(simu_tab, 15.87, axis=1)
    sv_minus2 = np.percentile(simu_tab, 2.275, axis=1)

    ax.plot(masses, null_hypothesis, label=r'$ H_{0} $ median', **kwargs)
    ax.plot(masses, sv_plus1, label=r'sv_plus1', **kwargs)
    ax.plot(masses, sv_plus2, label=r'sv_plus2', **kwargs)
    ax.plot(masses, sv_minus1, label=r'sv_minus1', **kwargs)
    ax.plot(masses, sv_minus2, label=r'sv_minus2', **kwargs)

    return ax
