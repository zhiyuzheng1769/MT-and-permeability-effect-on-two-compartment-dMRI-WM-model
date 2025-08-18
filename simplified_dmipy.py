'''Simiplified implementation of an isotropic gaussian diffusion model for modelling the extra cellular compartment in a two-compartment model. Adapted from dmipy (https://github.com/AthenaEPI/dmipy)'''

import numpy as np

water_gyromagnetic_ratio=267.513e6   # rad/(sT)

class Ball:
    r""" The Ball model [1]_ - an isotropic Tensor with one diffusivity.

    Parameters
    ----------
    lambda_iso : float,
        isotropic diffusivity in m^2/s.

    """
    def __init__(self, lambda_iso=None):
        self.lambda_iso = lambda_iso

    def __call__(self, acquisition_scheme, **kwargs):
        r'''
        Estimates the signal attenuation.

        Parameters
        ----------
        acquisition_scheme : AcquisitionScheme instance,
            An acquisition scheme that has been instantiated.
        Returns
        -------
        attenuation : float or array, shape(N),
            signal attenuation
        '''
        bvals = acquisition_scheme.bvalues
        lambda_iso = kwargs.get('lambda_iso', self.lambda_iso)
        E_ball = np.exp(-bvals * lambda_iso)
        return E_ball

class AcquisitionScheme:
    """
    Class that contains b values needed to simulate and
    fit data using microstructure models.
    """

    def __init__(self, bvalues):
        self.bvalues = bvalues.astype(float)

def acquisition_scheme_from_gradients(
        gradient_strengths, delta, Delta):
    r"""
    Creates an acquisition scheme object from gradient strengths, gradient
    directions pulse duration $\delta$ and pulse separation time $\Delta$.

    Parameters
    ----------
    gradient_strengths: 1D numpy array of shape (Ndata)
        gradient strength of the acquisition in T/m.
        e.g., a gradient strength of 300 mT/m must be entered as 300 / 1e3 T/m
    delta: float or 1D numpy array of shape (Ndata)
        if float, pulse duration of every measurements in seconds.
        if array, potentially varying pulse duration per measurement.
    Delta: float or 1D numpy array of shape (Ndata)
        if float, pulse separation time of every measurements in seconds.
        if array, potentially varying pulse separation time per measurement.

    Returns
    -------
    DmipyAcquisitionScheme: acquisition scheme object
        contains b values of the acquisition scheme to be used in any
        microstructure model.
    """
    bvalues = b_from_g(gradient_strengths, delta, Delta)
    return AcquisitionScheme(bvalues)

def b_from_g(
    g, delta, Delta,
    gyromagnetic_ratio=water_gyromagnetic_ratio
):
    """Compute b-value from gradient strength. Units are standard units."""
    tau = Delta - delta / 3
    return (g * gyromagnetic_ratio * delta) ** 2 * tau