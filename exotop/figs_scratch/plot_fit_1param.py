import numpy as np
import matplotlib.pyplot as plt
from useful_and_bespoke import colorize, colourbar, colourised_legend
from scipy import odr, stats
import pandas as pd
from matplotlib import rc
from matplotlib.pyplot import rcParams
from datetime import date
import matplotlib.lines as mlines

rc('text', usetex=True)
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'CMU Serif'
today = date.today().strftime("%b-%d-%Y")
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
labelsize = 22
legsize = 20
ticksize = 16

df = pd.read_csv('/home/claire/Works/aspect/runs/figs/aspect-output.csv')
print(df.head())
Ra = df.Ra_i_eff.to_numpy()
Ra1 = df.Ra_1.to_numpy()
h = df.h_rms.to_numpy()
eta = df.delta_eta.to_numpy()
ln_eta = np.log(eta)
eta_1 = np.exp(-ln_eta * 1)
eta_i = np.exp(-ln_eta * df.T_i.to_numpy())
Ra_i = Ra1 * eta_1 / eta_i
T_i = df.T_i.to_numpy()
alpha_m, dT_m, d_m = 2e-5, 1600, 2800e3
dim_factor = alpha_m * dT_m * d_m  # for testing

def printarr(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))


def func_lin(beta, u):
    a, b = beta
    return a*u + b


def func_lin_jac(u, **kwargs):
    # return jacobian for above function
    np.warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)
    try:
        return np.array([u, 1])
    except np.VisibleDeprecationWarning:
        return np.array([u, 1], dtype=object)  # this might not give the right answer


def fit_ODR(y, x1, x2=None, beta0=None, func=None, err_x=None, err_y=None, ci=0.95, maxit=1000):
    if x2 is not None:
        x = np.row_stack((x1, x2))  # odr doesn't seem to work with column_stack
    else:
        x = x1
    if beta0 is None:
        beta0 = [0.1, 1]
        if x2 is not None:
            beta0.append(1)
    data = odr.RealData(x, y, sx=err_x, sy=err_y)
    model = odr.Model(func)
    odrfit = odr.ODR(data, model, beta0, maxit=maxit)
    output = odrfit.run()

    # confidence intervals
    df_e = len(x1) - len(output.beta)  # degrees of freedom, error
    conf = []
    t_df = stats.t.ppf(ci, df_e)  # 0.975
    for i in range(len(output.beta)):
        conf.append([output.beta[i] - t_df * output.sd_beta[i],
                     output.beta[i] + t_df * output.sd_beta[i]])

    # chi sqr
    expected = func(output.beta, x)
    chisqr = output.sum_square  # np.sum(((y - expected) ** 2) / expected)
    chisqr_nu = output.res_var

    # covariance matrix
    cov_beta = output.cov_beta * output.res_var

    print('       -> ODR RESULTS')
    print('         -> reason for halting:', output.stopreason)
    for ii, val in enumerate(output.beta):
        print('         -> beta', ii, ':', val, '+/-', output.sd_beta[ii], '    CI:', conf[ii][0], conf[ii][1])
    print('         -> sum of squares error:', chisqr)
    print('         -> reduced chi sqr:', chisqr_nu)
    print('         -> covariance matrix:\n')
    printarr(cov_beta)
    # print('         -> eps', output.eps, len(output.eps), len(y))
    print('\n')
    return output


def fit_variance(cov_beta=None, jac_fn=None, xn=None, **kwargs):
    # variance of f(x1, x2) from error on fitted parameters, I think this won't work unless x1n and x2n are scalars

    jac = jac_fn(u=xn, **kwargs)
    # print('jac * cov =\n')
    # try:
    #     printarr(np.matmul(jac, cov_beta))
    # except TypeError:
    #     print(np.matmul(jac, cov_beta))

    var_y = np.matmul(np.matmul(jac, cov_beta), jac.T)
    # print('var_y', var_y)

    return var_y


def draw_random_yhat(u=None, beta=None, cov_beta=None, func=None):
    # beta_rand = [np.random.multivariate_normal(be, co) for (be, co) in zip(beta, cov_beta)]
    beta_rand = np.random.multivariate_normal(beta, cov_beta)
    yhat = func(beta_rand, u)
    return yhat


""" setup """
fig, ax = plt.subplots(figsize=(5,5))
beta0 = [0.01, -0.1]
x1 = np.log10(Ra)
y = np.log10(h)

# try fully linear model with interaction term
output = fit_ODR(y=y, x1=x1, beta0=beta0, func=func_lin, maxit=100000)
beta_h = output.beta
cov_beta = output.cov_beta * output.res_var

xlim = (6, 7.5)
x1_hat = np.linspace(xlim[0], xlim[-1])
h_hat = func_lin(beta_h, x1_hat)

c = colorize(ln_eta, cmap='rainbow', vmin=13, vmax=22)[0]
lneta_unique, idx = np.unique(ln_eta, return_index=True)
for n, b in enumerate(lneta_unique):
    mask = np.where(ln_eta == b)
    # plot original
    ax.plot(x1[mask], y[mask], c=c[idx[n]], lw=0, marker='o')

    # plot fit
    ax.plot(x1_hat, h_hat, c=c[idx[n]], lw=0.7)

    # error
    var_h = np.zeros_like(h_hat)
    for ii, x1n in enumerate(x1_hat):
        var_h[ii] = fit_variance(cov_beta=cov_beta, jac_fn=func_lin_jac, xn=x1n)
    SE_y = np.sqrt(var_h)
    yn_upper = h_hat + SE_y
    yn_lower = h_hat - SE_y
    ax.fill_between(x1_hat, yn_lower, yn_upper, fc=c[idx[n]], alpha=0.15)

    for nn in range(100):
        yhn = draw_random_yhat(u=x1_hat, beta=beta_h, cov_beta=cov_beta, func=func_lin)
        # ax.plot(x1_hat, yhn, c='k', alpha=0.1, lw=0.5)


handles = []
lw = [0, 0.7]
label = ['Data', 'Model']
ch = ('k', 'k')
marker = ['o', None]
for jj in range(2):
    handles.append(mlines.Line2D([], [], color=ch[jj], lw=lw[jj], label=label[jj], marker=marker[jj]))
leg = ax.legend(handles=handles, frameon=False, fontsize=legsize)
ax.add_artist(leg)

# colourbar(vector=lneta_unique, ax=ax, label='$b = \ln(\Delta \eta)$', labelsize=labelsize, labelpad=30)
ax = colourised_legend(clist=c, ax=ax, cleglabels=[int(e) for e in lneta_unique],
                       title='$\ln(\Delta \eta)$', legsize=ticksize, titlesize=legsize, markersize=5)
ax.set_xlim(xlim)
# ax.set_ylim(-2.25, -1.95)

ax.set_xlabel(r'log(Ra$_{i, eff}$)', fontsize=labelsize)
ax.set_ylabel(r'log(h$^\prime_{\rm rms}$)', fontsize=labelsize)
ax.tick_params(axis='x', labelsize=ticksize)
ax.tick_params(axis='y', labelsize=ticksize)

fig.savefig(fig_path + 'hfit_1_' + today + '.png', bbox_inches='tight')
plt.show()
