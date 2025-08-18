# Import required libraries for diffusion modeling, optimization, and data handling
import numpy as np  
from scipy.optimize import minimize  # Optimization tools
import pandas as pd  
from scipy.interpolate import CubicSpline  # Spline interpolation
import sys
sys.path.append("..")
from simplified_dmipy import Ball, acquisition_scheme_from_gradients

# Function to construct the AxCaliber forward model
def forge_axcaliber():
    # Recreate the sequences used in the original AxCaliber paper: https://doi.org/10.1002/mrm.21577
    gradient_strengths = np.squeeze(np.tile(np.arange(0, 1.21, 0.08), (1, 8)))  # Gradient strengths (T/m)

    # Define timing parameters for the acquisition scheme
    delta = 0.0025  # Gradient duration (s)
    Delta = np.tile(np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08]), (16, 1)).transpose().flatten()  # Diffusion times (s)
    
    # Create the acquisition scheme
    acq_scheme = acquisition_scheme_from_gradients(gradient_strengths, delta, Delta)
    # Define a Gaussian "Ball" model (for isotropic diffusion) that represents extraaxonal compartment
    ball = Ball()

    # Construct the interpolated DW attenuation dictionary for different radii
    # Initialize a list to store diffusion weighted signal attenuations at different radius values
    ys = [np.ones([128])]  # Start with a uniform array of ones for r=0

    #Loop through all radius values and read the signal and calculate the intraaxonal DW signal attenuation and append the values to the dictionary
    for i in np.arange(0.30, 5.51, 0.10):  # Radii in micrometers
        # Read simulated signal data for a specific radius
        simdata = pd.read_csv("./perm_results/signal_MT_0_sus_0_perm_0.000_rmean_" + '{0:.2f}'.format(i) + "_density_0.65.csv", header=None)
        
        # calculate the intraaxonal DW signal attenuation by dividing the signal magnitude by the b0 signal magnitude
        y = simdata[0].to_numpy() / simdata[0][0]
        ys.append(y)  # Append normalized signal to the list

    # Convert the list of signals into a NumPy array for fitting
    np_ys = np.array(ys)
    
    # Create an array of radii (including 0 at the start)
    rs = np.arange(0.30, 5.51, 0.10)
    rs = np.insert(rs, 0, 0)  # Insert radius = 0 at the beginning

    # Create a cubic spline interpolation for the signal as a function of radius
    spline_cyl = CubicSpline(rs, np_ys)

    # Define the AxCaliber two-compartment model
    def axcaliber(fr, Dh, r):
        # Combine the contributions of the ball for extraaxonal signal and the spline interpolated dictionary for intraaxonal signal
        return (1 - fr) * ball(acq_scheme, lambda_iso=Dh) + fr * spline_cyl(r)

    return axcaliber

# Basic function to fit the AxCaliber forward model to data, y is measured/simulated data to be fitted, fr, Dh, r are initial guesses for the fitting algorithm, forward_model is the model we will use to fit the data
def fit_dict_axcaliber(y, fr, Dh, r, forward_model):
    # Define the cost function for optimization
    def cost_func(params):
        f, Dh, r = params  # Unpack parameters, f=volume fraction, Dh=extraaxonal diffusivity (um^2/ms), r=(mean) cylinder radius (um)
        # Compute the cost as the sum of squared error between model prediction and data
        result = np.sum(((forward_model(f, Dh * 1e-9, r) - y)) ** 2) # Note the Dh is defined in um^2/ms, but dmipy uses m^2/s
        return result

    # Initial guess for parameters: volume fraction (fr), diffusivity (Dh), and radius (r)
    initial_params = [fr, Dh, r]

    # Use the Nelder-Mead method to minimize the cost function
    popt = minimize(fun=cost_func, x0=initial_params, method='Nelder-Mead',
                    bounds=[(0, 1), (0., 4), (0., 5.5)])  # Parameter bounds for f, Dh, r respectively

    return popt  # Return optimized parameters

# The overall fitting function that calls the fit_dict_axcaliber() above with a grid of different initial guesses and return the fitting result that has minimum error
def fit_params(forward_model, y):
    # Define optimization methods to try
    optimisers = ['Nelder-Mead']

    # Lists to store results of optimization
    nm = []  # Store parameter sets
    nm_cost = []  # Store corresponding cost values

    # Loop over initial guesses for each parameter
    for fr in np.arange(0.05, 1, 0.45):  # Volume fraction initial guesses
        for r in np.arange(0.5, 5, 2):  # Radius initial guesses
            for d in np.arange(0.2, 2, 0.8):  # Diffusivity initial guesses
                # Perform optimization with current initial guess
                nm_fit = fit_dict_axcaliber(y, fr, d, r, forward_model)
                nm.append(nm_fit.x)  # Append the optimized parameters
                nm_cost.append(nm_fit.fun)  # Append the corresponding cost

                lis = [nm]  # Wrap the parameters in a list

    # Find the minimum cost and the best parameters
    min_costs = [np.min(nm_cost)]  # Minimum cost value
    opts = [lis[0][np.argmin(nm_cost)]]  # Best parameter set
    
    # Create dictionaries to store results for each optimizer
    min_cost_dict = dict(zip(optimisers, min_costs))
    opts_dict = dict(zip(optimisers, opts))

    return min_cost_dict, opts_dict  # Return results
