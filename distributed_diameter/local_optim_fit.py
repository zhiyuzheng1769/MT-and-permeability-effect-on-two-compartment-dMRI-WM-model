from dmipy.signal_models import gaussian_models, cylinder_models
from dmipy.core.acquisition_scheme import acquisition_scheme_from_qvalues
import numpy as np
from scipy.optimize import curve_fit, minimize
import pandas as pd
from scipy.stats import gamma
from dmipy.core.acquisition_scheme import acquisition_scheme_from_gradient_strengths
from scipy.interpolate import CubicSpline

def reorder_df(dataframe):
    new_df = dataframe.copy(deep=True)
    for i in range(dataframe.shape[0]):
        if i % 16 == 1:
            new_df.loc[i + 12] = dataframe.loc[i]
            new_df.loc[i] = dataframe.loc[i + 12]
        elif i % 16 > 1 and i % 16 <= 3:
            new_df.loc[i + 12] = dataframe.loc[i]
            new_df.loc[i] = dataframe.loc[i + 2]
        elif i % 16 > 3 and i % 16 < 11:
            new_df.loc[i] = dataframe.loc[i + 2]
        elif i % 16 >= 11 and i % 16 < 13:
            new_df.loc[i] = dataframe.loc[i + 3]
    return new_df

# def forge_axcaliber(method="curve_fit", y=None):
#     Dr = 3e-9
#     gradient_strengths = np.squeeze(np.tile(np.arange(0, 1.21, 0.08), (1,8)) )  # T/m
#     gradient_directions = np.tile(np.array([0,1,0]), (128, 1))# np.asarray([[0, 0, 1],[0, 0, 1],[0, 0, 1]])  # Perpendicular to the cylinder
#     delta = 0.0025  # s
#     Delta = np.tile(np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08]), (16,1)).transpose().flatten()  # s
#     acq_scheme = acquisition_scheme_from_gradient_strengths(gradient_strengths, gradient_directions, delta, Delta)
#     ball = gaussian_models.G1Ball()
#     ys = [np.ones([128])]
#     for i in np.arange(0.30, 5.51, 0.10):
#         # a = '{0:.2f}'.format(i)
#         simdata = pd.read_csv("./perm/signal_MT_0_sus_0_perm_0.000_rmean_" + '{0:.2f}'.format(i) + "_density_0.65.csv",header=None)

#         y = simdata[0].to_numpy()/simdata[0][0]
#         ys.append(y)
#     np_ys = np.array(ys)
#     rs = np.arange(0.30, 5.51, 0.10)
#     rs = np.insert(rs, 0, 0)
#     spline_cyl = CubicSpline(rs, np_ys)
#     def axcaliber(fr, Dh, r):
#         return (1 - fr)*ball(acq_scheme, lambda_iso=Dh) + fr*spline_cyl(r)
#     return axcaliber
    
def forge_axcaliber_250k(method="curve_fit", y=None):
    Dr = 3e-9
    gradient_strengths = np.squeeze(np.tile(np.arange(0, 1.21, 0.08), (1,8)) )  # T/m
    gradient_directions = np.tile(np.array([0,1,0]), (128, 1))# np.asarray([[0, 0, 1],[0, 0, 1],[0, 0, 1]])  # Perpendicular to the cylinder
    delta = 0.0025  # s
    Delta = np.tile(np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08]), (16,1)).transpose().flatten()  # s
    acq_scheme = acquisition_scheme_from_gradient_strengths(gradient_strengths, gradient_directions, delta, Delta)
    ball = gaussian_models.G1Ball()
    ys = [np.ones([128])]
    for i in np.arange(0.30, 5.51, 0.10):
        # a = '{0:.2f}'.format(i)
        simdata = pd.read_csv("../same_packing_test/same_packing_scaled_bound/perm/signal_MT_0_sus_0_perm_0.000_rmean_" + '{0:.2f}'.format(i) + "_density_0.65.csv",header=None)

        y = simdata[0].to_numpy()/simdata[0][0]
        ys.append(y)
    np_ys = np.array(ys)
    rs = np.arange(0.30, 5.51, 0.10)
    rs = np.insert(rs, 0, 0)
    spline_cyl = CubicSpline(rs, np_ys)
    def axcaliber(fr, Dh, r):
        return (1 - fr)*ball(acq_scheme, lambda_iso=Dh) + fr*spline_cyl(r)
    return axcaliber

def fit_dict_axcaliber(y, fr, Dh, r, forward_model, method="curve_fit", cost="LSE"):
    
    def forward_curve_fit(x, *params):
        return forward_model(*params)
    def cost_func(params):
        f, Dh, r = params
        result = np.sum(((forward_model(f, Dh * 1e-9, r) - y)) ** 2)
        return result
    def cost_func_div(params):
        f, Dh, r = params
        result = np.sum(((forward_model(f, Dh * 1e-9, r) - y)/y) ** 2)
        return result
    if method == "curve_fit":
        initial_params = [fr, Dh*1e-9, r]
        popt, pcov = curve_fit(f=forward_curve_fit, xdata=0, ydata=y, p0=initial_params, maxfev=5000, bounds=([0, 0, 0.], [1, 4e-9, 5.5]))
        popt[1] = popt[1]*1e9
    else:
        initial_params = [fr, Dh, r]
        if cost == "LSE":
            popt = minimize(fun=cost_func, x0=initial_params, method=method, bounds=[(0, 1), (0., 4), (0., 5.5)])
        elif cost == "div":
            popt = minimize(fun=cost_func_div, x0=initial_params, method=method, bounds=[(0, 1), (0., 4), (0., 5.5)])
    return popt

def local_optim_test_indata(forward_model, y, cost="LSE"):    
    optimisers = ['curve_fit','Nelder-Mead','Powell','L-BFGS-B','TNC','COBYLA','SLSQP','trust-constr']
    cf = []
    nm = []
    pl = []
    lb = []
    tnc = []
    cobyla = []
    slsqp = []
    tc = []
    cf_cost = []
    nm_cost = []
    pl_cost = []
    lb_cost = []
    tnc_cost = []
    cobyla_cost = []
    slsqp_cost = []
    tc_cost = []
    # forward_model = forge_axcaliber()
    for fr in np.arange(0.05, 1, 0.45):
        for r in np.arange(0.5, 5, 2):
            for d in np.arange(0.2, 2, 0.8):
            
                cf_fit = fit_dict_axcaliber(y, fr, d, r, forward_model, method="curve_fit", cost=cost)
                nm_fit = fit_dict_axcaliber(y, fr, d, r, forward_model, method="Nelder-Mead", cost=cost)
                pl_fit = fit_dict_axcaliber(y, fr, d, r, forward_model, method="Powell", cost=cost)
                lb_fit = fit_dict_axcaliber(y, fr, d, r, forward_model, method="L-BFGS-B", cost=cost)
                tnc_fit = fit_dict_axcaliber(y, fr, d, r,forward_model, method="TNC", cost=cost)
                cobyla_fit = fit_dict_axcaliber(y, fr, d, r, forward_model, method="COBYLA", cost=cost)
                slsqp_fit = fit_dict_axcaliber(y, fr, d, r, forward_model, method="SLSQP", cost=cost)
                tc_fit = fit_dict_axcaliber(y, fr, d, r, forward_model, method="trust-constr", cost=cost)
                cf.append(cf_fit)
                cost_cf = np.sum((forward_model(*cf_fit) - y) ** 2)
                cf_cost.append(cost_cf)
                nm.append(nm_fit.x)
                pl.append(pl_fit.x)
                lb.append(lb_fit.x)
                tnc.append(tnc_fit.x)
                cobyla.append(cobyla_fit.x)
                slsqp.append(slsqp_fit.x)
                tc.append(tc_fit.x)
                nm_cost.append(nm_fit.fun)
                pl_cost.append(pl_fit.fun)
                lb_cost.append(lb_fit.fun)
                tnc_cost.append(tnc_fit.fun)
                cobyla_cost.append(cobyla_fit.fun)
                slsqp_cost.append(slsqp_fit.fun)
                tc_cost.append(tc_fit.fun)
                lis = [cf,nm,pl,lb,tnc,cobyla,slsqp,tc]
    
    min_costs = [np.min(cf_cost),np.min(nm_cost),np.min(pl_cost),np.min(lb_cost),np.min(tnc_cost),np.min(cobyla_cost),np.min(slsqp_cost),np.min(tc_cost)]
    opts = [lis[0][np.argmin(cf_cost)],lis[1][np.argmin(nm_cost)],lis[2][np.argmin(pl_cost)],lis[3][np.argmin(lb_cost)],lis[4][np.argmin(tnc_cost)],lis[5][np.argmin(cobyla_cost)],lis[6][np.argmin(slsqp_cost)],lis[7][np.argmin(tc_cost)]]
    min_cost_dict = dict(zip(optimisers, min_costs))
    opts_dict = dict(zip(optimisers, opts))
    return min_cost_dict, opts_dict

def forge_axcaliber_semifinal_dt(method="curve_fit", y=None):
    Dr = 3e-9
    gradient_strengths = np.squeeze(np.tile(np.arange(0, 1.21, 0.08), (1,8)) )  # T/m
    gradient_directions = np.tile(np.array([0,1,0]), (128, 1))# np.asarray([[0, 0, 1],[0, 0, 1],[0, 0, 1]])  # Perpendicular to the cylinder
    delta = 0.0025  # s
    Delta = np.tile(np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08]), (16,1)).transpose().flatten()  # s
    acq_scheme = acquisition_scheme_from_gradient_strengths(gradient_strengths, gradient_directions, delta, Delta)
    ball = gaussian_models.G1Ball()
    ys = [np.ones([128])]
    for i in np.arange(0.30, 5.51, 0.10):
        # a = '{0:.2f}'.format(i)
        simdata = pd.read_csv("./100ktest/signal_MT_0_sus_0_perm_0.000_rmean_" + '{0:.2f}'.format(i) + "_density_0.65.csv",header=None)

        y = simdata[0].to_numpy()/simdata[0][0]
        ys.append(y)
    np_ys = np.array(ys)
    rs = np.arange(0.30, 5.51, 0.10)
    rs = np.insert(rs, 0, 0)
    spline_cyl = CubicSpline(rs, np_ys)
    def axcaliber(fr, Dh, r, A):
        return (1 - fr)*np.exp(-acq_scheme.bvalues*(Dh + A*1e-12*(np.log(Delta/delta) + 1.5)/(Delta - delta/3))) + fr*spline_cyl(r)
    return axcaliber

def local_optim_test_indata_semifinal_dt(forward_model, y, cost="LSE"):    
    optimisers = ['curve_fit','Nelder-Mead','Powell','L-BFGS-B','TNC','COBYLA','SLSQP','trust-constr']
    cf = []
    nm = []
    pl = []
    lb = []
    tnc = []
    cobyla = []
    slsqp = []
    tc = []
    cf_cost = []
    nm_cost = []
    pl_cost = []
    lb_cost = []
    tnc_cost = []
    cobyla_cost = []
    slsqp_cost = []
    tc_cost = []
    # forward_model = forge_axcaliber()
    for fr in np.arange(0.05, 1, 0.45):
        for r in np.arange(0.5, 5, 2):
            for d in np.arange(0.2, 2, 0.8):
                for A in [-3, 0, 3]:
            
                    # cf_fit = fit_dict_axcaliber_dt(y, fr, d, r, A, forward_model, method="curve_fit")
                    nm_fit = fit_dict_axcaliber_dt(y, fr, d, r, A, forward_model, method="Nelder-Mead")
                    pl_fit = fit_dict_axcaliber_dt(y, fr, d, r, A, forward_model, method="Powell")
                    lb_fit = fit_dict_axcaliber_dt(y, fr, d, r, A, forward_model, method="L-BFGS-B")
                    tnc_fit = fit_dict_axcaliber_dt(y, fr, d, r, A,forward_model, method="TNC")
                    cobyla_fit = fit_dict_axcaliber_dt(y, fr, d, r, A, forward_model, method="COBYLA")
                    slsqp_fit = fit_dict_axcaliber_dt(y, fr, d, r, A, forward_model, method="SLSQP")
                    tc_fit = fit_dict_axcaliber_dt(y, fr, d, r, A, forward_model, method="trust-constr")
                    # cf.append(cf_fit)
                    # cost_cf = np.sum((forward_model(*cf_fit) - y) ** 2)
                    # cf_cost.append(cost_cf)
                    nm.append(nm_fit.x)
                    pl.append(pl_fit.x)
                    lb.append(lb_fit.x)
                    tnc.append(tnc_fit.x)
                    cobyla.append(cobyla_fit.x)
                    slsqp.append(slsqp_fit.x)
                    tc.append(tc_fit.x)
                    nm_cost.append(nm_fit.fun)
                    pl_cost.append(pl_fit.fun)
                    lb_cost.append(lb_fit.fun)
                    tnc_cost.append(tnc_fit.fun)
                    cobyla_cost.append(cobyla_fit.fun)
                    slsqp_cost.append(slsqp_fit.fun)
                    tc_cost.append(tc_fit.fun)
                    lis = [nm,pl,lb,tnc,cobyla,slsqp,tc]
    
    min_costs = [np.min(nm_cost),np.min(pl_cost),np.min(lb_cost),np.min(tnc_cost),np.min(cobyla_cost),np.min(slsqp_cost),np.min(tc_cost)]
    opts = [lis[0][np.argmin(nm_cost)],lis[1][np.argmin(pl_cost)],lis[2][np.argmin(lb_cost)],lis[3][np.argmin(tnc_cost)],lis[4][np.argmin(cobyla_cost)],lis[5][np.argmin(slsqp_cost)],lis[6][np.argmin(tc_cost)]]
    min_cost_dict = dict(zip(optimisers, min_costs))
    opts_dict = dict(zip(optimisers, opts))
    return min_cost_dict, opts_dict

def forge_axcaliber_kurt():
    Dr = 3e-9
    gradient_strengths = np.squeeze(np.tile(np.arange(0, 1.21, 0.08), (1,8)) )  # T/m
    gradient_directions = np.tile(np.array([0,1,0]), (128, 1))# np.asarray([[0, 0, 1],[0, 0, 1],[0, 0, 1]])  # Perpendicular to the cylinder
    delta = 0.0025  # s
    Delta = np.tile(np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08]), (16,1)).transpose().flatten()  # s
    acq_scheme = acquisition_scheme_from_gradient_strengths(gradient_strengths, gradient_directions, delta, Delta)
    ball = gaussian_models.G1Ball()
    ys = [np.ones([128])]
    for i in np.arange(0.30, 5.51, 0.10):
        # a = '{0:.2f}'.format(i)
        simdata = pd.read_csv("./100ktest/signal_MT_0_sus_0_perm_0.000_rmean_" + '{0:.2f}'.format(i) + "_density_0.65.csv",header=None)

        y = simdata[0].to_numpy()/simdata[0][0]
        ys.append(y)
    np_ys = np.array(ys)
    rs = np.arange(0.30, 5.51, 0.10)
    rs = np.insert(rs, 0, 0)
    spline_cyl = CubicSpline(rs, np_ys)
    def axcaliber(fr, Dh, r, K):
        return (1 - fr)*np.exp(-acq_scheme.bvalues*Dh + (acq_scheme.bvalues**2)*(Dh**2)*K/6) + fr*spline_cyl(r)
    return axcaliber

def fit_dict_axcaliber_kurt(y, fr, Dh, r, K, forward_model, method="curve_fit"):
    
    def forward_curve_fit(x, *params):
        # print(*params)
        return forward_model(*params)
    def cost_func(params):
        f, Dh, r, K = params
        # print(params)
        result = np.sum((forward_model(f, Dh * 1e-9, r, K) - y) ** 2)
        return result
    if method == "curve_fit":
        initial_params = [fr, Dh*1e-9, r, K]
        popt, pcov = curve_fit(f=forward_curve_fit, xdata=0, ydata=y, p0=initial_params, maxfev=5000, bounds=([0, 0, 0., 0], [1, 4e-9, 10., np.inf]))
        popt[1] = popt[1]*1e9
    else:
        initial_params = [fr, Dh, r, K]
        popt = minimize(fun=cost_func, x0=initial_params, method=method, bounds=[(0, 1), (0., 4), (0., 5.5), (0, np.inf)])
    return popt

def local_optim_test_indata_kurt(forward_model, y):    
    optimisers = ['curve_fit','Nelder-Mead','Powell','L-BFGS-B','TNC','COBYLA','SLSQP','trust-constr']
    cf = []
    nm = []
    pl = []
    lb = []
    tnc = []
    cobyla = []
    slsqp = []
    tc = []
    cf_cost = []
    nm_cost = []
    pl_cost = []
    lb_cost = []
    tnc_cost = []
    cobyla_cost = []
    slsqp_cost = []
    tc_cost = []
    K = 0.1
    # forward_model = forge_axcaliber()
    for fr in np.arange(0.05, 1, 0.45):
        for r in np.arange(0.5, 5, 2):
            for d in np.arange(0.2, 2, 0.8):
                for K in np.arange(0.1, 1, 0.4):
                    cf_fit = [] # fit_dict_axcaliber_kurt(y, fr, d, r, K, forward_model, method="curve_fit")
                    nm_fit = fit_dict_axcaliber_kurt(y, fr, d, r, K, forward_model, method="Nelder-Mead")
                    pl_fit = fit_dict_axcaliber_kurt(y, fr, d, r, K, forward_model, method="Powell")
                    lb_fit = fit_dict_axcaliber_kurt(y, fr, d, r, K, forward_model, method="L-BFGS-B")
                    tnc_fit = fit_dict_axcaliber_kurt(y, fr, d, r, K,forward_model, method="TNC")
                    cobyla_fit = fit_dict_axcaliber_kurt(y, fr, d, r, K, forward_model, method="COBYLA")
                    slsqp_fit = fit_dict_axcaliber_kurt(y, fr, d, r, K, forward_model, method="SLSQP")
                    tc_fit = fit_dict_axcaliber_kurt(y, fr, d, r, K, forward_model, method="trust-constr")
                    cf.append(cf_fit)
                    # cost_cf = np.sum((forward_model(*cf_fit) - y) ** 2)
                    cf_cost.append(0)
                    nm.append(nm_fit.x)
                    pl.append(pl_fit.x)
                    lb.append(lb_fit.x)
                    tnc.append(tnc_fit.x)
                    cobyla.append(cobyla_fit.x)
                    slsqp.append(slsqp_fit.x)
                    tc.append(tc_fit.x)
                    nm_cost.append(nm_fit.fun)
                    pl_cost.append(pl_fit.fun)
                    lb_cost.append(lb_fit.fun)
                    tnc_cost.append(tnc_fit.fun)
                    cobyla_cost.append(cobyla_fit.fun)
                    slsqp_cost.append(slsqp_fit.fun)
                    tc_cost.append(tc_fit.fun)
                    lis = [cf,nm,pl,lb,tnc,cobyla,slsqp,tc]
    
    min_costs = [np.min(cf_cost),np.min(nm_cost),np.min(pl_cost),np.min(lb_cost),np.min(tnc_cost),np.min(cobyla_cost),np.min(slsqp_cost),np.min(tc_cost)]
    opts = [lis[0][np.argmin(cf_cost)],lis[1][np.argmin(nm_cost)],lis[2][np.argmin(pl_cost)],lis[3][np.argmin(lb_cost)],lis[4][np.argmin(tnc_cost)],lis[5][np.argmin(cobyla_cost)],lis[6][np.argmin(slsqp_cost)],lis[7][np.argmin(tc_cost)]]
    min_cost_dict = dict(zip(optimisers, min_costs))
    opts_dict = dict(zip(optimisers, opts))
    return min_cost_dict, opts_dict