import matplotlib.pyplot as plt
import numpy as np
import torch
from IPython.display import display, clear_output 


def bishop_funcs(x, heaviside_x0 = 0):
    """ 
    Compute functions used in Bishop chapter 5 to illustrate network function estimation capabilities.
        - x : an array of values for which the functions are computed
        - heaviside_x0 : center position for the heaviside function
    """
    
    y_a = x**2
    y_b = np.sin(x)
    y_c = np.abs(x)
    y_d = np.heaviside(x, heaviside_x0)
    
    y = {"a": y_a, "b": y_b, "c": y_c, "d": y_d, "x":x}
    
    return y


def plot_bishop_func(y, y_n = None):
    
    fig = plt.figure(figsize=(8,8), constrained_layout=True) # Initiate figure with constrained layout
    gs = fig.add_gridspec(2, 2) # Add 2x2 grid

    ax1 = fig.add_subplot(gs[0, 0]) 
    ax1.set_title('(a)')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.grid()
    im1 = ax1.plot(y["x"], y["a"], '.')

    ax2 = fig.add_subplot(gs[0, 1]) 
    ax2.set_title('(b)')
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.grid()
    im2 = ax2.plot(y["x"], y["b"], '.')

    ax3 = fig.add_subplot(gs[1, 0]) 
    ax3.set_title('(c)')
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.grid()
    im3 = ax3.plot(y["x"], y["c"], '.')

    ax4 = fig.add_subplot(gs[1, 1]) 
    ax4.set_title('(d)')
    ax4.set_xlabel("x")
    ax4.set_ylabel("y")
    ax4.grid()
    im4 = ax4.plot(y["x"], y["d"], '.')
    
    if y_n is dict:
        ax1.plot(y_n["x"], y_n["a"], '-')
        ax2.plot(y_n["x"], y_n["b"], '-')
        ax3.plot(y_n["x"], y_n["c"], '-')
        ax4.plot(y_n["x"], y_n["d"], '-')


def plot_bishop_init():
    fig = plt.figure(figsize=(8,8), constrained_layout=True) # Initiate figure with constrained layout
    gs = fig.add_gridspec(3, 2) # Add 2x2 grid
    ax0 = fig.add_subplot(gs[0, :]) 
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[2, 1])
    
    return fig, gs, ax0, ax1, ax2, ax3, ax4


def plot_bishop_update(y_out, y_bishop, x_in, epoch, E_collect, params, fig, gs, ax0, ax1, ax2, ax3, ax4, include_weights = False):

    y_n = {"x": x_in.detach().numpy()}
    
    if include_weights == True:
        params_a, params_b, params_c, params_d = torch.chunk(params, 4, dim=0)

    y_na, y_nb, y_nc, y_nd = torch.chunk(y_out, 4, dim=-1)
    
    y_n["a"] = y_na.detach().numpy()
    y_n["b"] = y_nb.detach().numpy()
    y_n["c"] = y_nc.detach().numpy()
    y_n["d"] = y_nd.detach().numpy()
    
    ax0.clear()
    ax0.set_title('Error')
    ax0.set_xlabel("Epoch")
    ax0.set_ylabel("Sum of MSE")
    ax0.grid()
    im0 = ax0.semilogy(np.arange(1,epoch+1), E_collect, '-', color="C0")

    ax1.clear() 
    ax1.set_title('(a)')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.grid()
    if include_weights == True:
        ax1.plot(y_n["x"], params_a.detach().numpy()[:,0], '--', color="C2", linewidth = 0.5)
        ax1.plot(y_n["x"], params_a.detach().numpy()[:,1], '--', color="C3", linewidth = 0.5)
        ax1.plot(y_n["x"], params_a.detach().numpy()[:,2], '--', color="C4", linewidth = 0.5)
    ax1.plot(y_n["x"], y_n["a"], '-', color="C1", linewidth = 1)
    im1 = ax1.plot(y_bishop["x"], y_bishop["a"], '.', color="C0", markersize = 3.0)
 
    ax2.clear() 
    ax2.set_title('(b)')
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.grid()
    if include_weights == True:
        ax2.plot(y_n["x"], params_b.detach().numpy()[:,0], '--', color="C2", linewidth = 0.5)
        ax2.plot(y_n["x"], params_b.detach().numpy()[:,1], '--', color="C3", linewidth = 0.5)
        ax2.plot(y_n["x"], params_b.detach().numpy()[:,2], '--', color="C4", linewidth = 0.5)
    ax2.plot(y_n["x"], y_n["b"], '-', color="C1", linewidth = 1)
    im2 = ax2.plot(y_bishop["x"], y_bishop["b"], '.', color="C0", markersize = 3.0)
 
    ax3.clear() 
    ax3.set_title('(c)')
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.grid()
    if include_weights == True:
        ax3.plot(y_n["x"], params_c.detach().numpy()[:,0], '--', color="C2", linewidth = 0.5)
        ax3.plot(y_n["x"], params_c.detach().numpy()[:,1], '--', color="C3", linewidth = 0.5)
        ax3.plot(y_n["x"], params_c.detach().numpy()[:,2], '--', color="C4", linewidth = 0.5)
    ax3.plot(y_n["x"], y_n["c"], '-', color="C1", linewidth = 1)
    im3 = ax3.plot(y_bishop["x"], y_bishop["c"], '.', color="C0", markersize = 3.0)

    ax4.clear() 
    ax4.set_title('(d)')
    ax4.set_xlabel("x")
    ax4.set_ylabel("y")
    ax4.grid()
    if include_weights == True:
        ax4.plot(y_n["x"], params_d.detach().numpy()[:,0], '--', color="C2", linewidth = 0.5)
        ax4.plot(y_n["x"], params_d.detach().numpy()[:,1], '--', color="C3", linewidth = 0.5)
        ax4.plot(y_n["x"], params_d.detach().numpy()[:,2], '--', color="C4", linewidth = 0.5)
    ax4.plot(y_n["x"], y_n["d"], '-', color="C1", linewidth = 1)
    im4 = ax4.plot(y_bishop["x"], y_bishop["d"], '.', color="C0", markersize = 3.0)    
    
    # To display result at each iteration
    display(fig)
    clear_output(wait=True)


def printProgressBar (iteration, total, *args, subject='', prefix = '', suffix = '', decimals = 1, length = 10, fill = 'O'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    
    if args:
        print('\r%s |%s| %s%% %s %s. Counter: %s/%s, Running error magnitude: %.1f' % (prefix, bar, percent, suffix, subject, iteration, total, args[0]), end = '\r')
    else:
        print('\r%s |%s| %s%% %s %s. Counter: %s/%s' % (prefix, bar, percent, suffix, subject, iteration, total), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()