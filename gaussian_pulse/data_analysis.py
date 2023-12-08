import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import cmath
import lmfit
from joblib import Parallel, delayed

from numpy.fft import fft, fftfreq, ifft

def band_filter(x):
    ar = []
    for i in x:
        if i >= 0.3 and i <=0.5:
            ar.append(1)
        else:
            ar.append(0)
    return np.array(ar)

def filter_data(time, data):
    filtered_data = []
    w0 = 2.35
    for j in range(data.shape[1]):
        signal = data[:,j]*np.exp(1j*w0*time)
        signal_fft = fft(signal)
        x_freq = fftfreq(len(time), d=1)
        filter_signal = signal_fft * band_filter(x_freq)
        filtered_data.append(ifft(filter_signal)*np.exp(-1j*w0*time))
    return np.transpose(np.array(filtered_data))

def convert_to_complex(df, i, j):
    return complex((df[j][i][:-2] + 'j').replace("*^","e"))

def extract_data(mat: str, kyvalues: list):
    theory_df00 = pd.read_csv('./'+mat+'/kx_0_ky_0/datafile.csv', header=None)
    theory_df01 = pd.read_csv('./'+mat+'/kx_0_ky_1/datafile.csv', header=None)
    theory_df02 = pd.read_csv('./'+mat+'/kx_0_ky_2/datafile.csv', header=None)
    theory_df03 = pd.read_csv('./'+mat+'/kx_0_ky_3/datafile.csv', header=None)
    theory_df04 = pd.read_csv('./'+mat+'/kx_0_ky_4/datafile.csv', header=None)
    theory_df05 = pd.read_csv('./'+mat+'/kx_0_ky_5/datafile.csv', header=None)

    data_row, data_column = theory_df00.shape
    time_domain = theory_df00.iloc[:,0]

    theory_data = np.zeros((data_row, data_column), dtype=complex)
    for i in range(data_row):
        theory_data[i][0] = theory_df00[0][i]
        for j in range(1, data_column):
            # theory_data[i][j] = (convert_to_complex(theory_df00, i, j))
            theory_data[i][j] = 0
            for index in kyvalues:
                theory_data[i][j] += convert_to_complex(eval("theory_df0" + str(index)), i, j)
    w0 = 2.35
    e_data = 1j*w0*theory_data[:, 1:]
    e_data = filter_data(time_domain, e_data)

    abs_data = np.zeros((data_row, data_column-1))
    arg_data = np.zeros((data_row, data_column-1))
    for i in range(data_row):
        for j in range(data_column-1):
            (abs_data[i,j], arg_data[i,j]) = cmath.polar(e_data[i,j])

    summed_abs_data = np.zeros(data_column-1)
    for i in range(data_column-1):
        summed_abs_data[i] = np.sum(abs_data[:,i])
    return time_domain.values, abs_data, arg_data, summed_abs_data

def max_E_index(mat: str, kyvalues: list, i):
    _, abs_data, _, _ = extract_data(mat, kyvalues)
    return max(range(len(abs_data[:,i])), key = abs_data[:,i].__getitem__)

def real_time_plot(mat: str, kyvalues: list, delay_index: int):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    fig.set_size_inches(10, 8)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    time_array, abs_data, _, _ = extract_data(mat, kyvalues)
    

    # plt.plot([time_data[i][0] for i in range(len(time_data.axes[1]))], [E_data1_45[i][delay_index] for i in range(len(E_data1_45.axes[1]))], '-o')
    ax.plot(time_array - time_array[max_E_index(mat, kyvalues, delay_index)], abs_data[:,delay_index]) 
    #Maxinmum value of E shifted to t=0
    # plt.legend(['Experimental data', 'Numerical calculation'],fontsize=18)
    plt.xlabel("Time (fs)", fontsize=25)
    plt.ylabel("E-field intensity", fontsize=25)
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    # ax.set_yticks([0, 1.4*(10**6)])
    plt.xlim([-400,400])
    # plt.show()
    # plt.savefig('./graphics/e-field.svg')
    return ax

def time_delay_plot(mat: str, kyvalues: list):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    fig.set_size_inches(10, 8)
    # ax.spines['left'].set_position('center')
    # # ax.spines['bottom'].set_position('zero')
    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')
    # ax.xaxis.set_ticks_position('bottom')
    # ax.yaxis.set_ticks_position('left')

    _, _, _, summed_abs_data = extract_data(mat, kyvalues)
    ax.plot(np.arange(10, -10.2, -0.2), summed_abs_data,"-*")
    # plt.plot(np.arange(30,-30.2,-0.2), correction_factor_02*summed_E_data_02)
    # plt.xlabel("Time Delay (fs)", fontsize=25)
    # plt.ylabel("Summed E-field ", fontsize=25)
    # ax.tick_params(axis='x', labelsize=25)
    # ax.tick_params(axis='y', labelsize=25)
    # plt.legend(['Numerical Data'],fontsize=18)
    # plt.xlim([-10.2,10])
    return ax
    # plt.show()

def phase_fun(x, phi_0, a, b, c, d):
        return phi_0 + a*x + b*(x)**2 + c*(x)**3 + d*(x)**4

def phase_fitting(mat: str, kyvalues: list, delay_index: int, low_lim: float, high_lim: float):
    model = lmfit.Model(phase_fun)
    time_array, abs_data, arg_data, _ = extract_data(mat, kyvalues)

    max_field_index = max_E_index(mat, kyvalues, delay_index)
    x = time_array - time_array[max_field_index]
    delaytime_range = np.logical_and(x>=low_lim, x<=high_lim)
    x1 = x[delaytime_range]

    y = np.unwrap(arg_data[:,delay_index]) - np.unwrap(arg_data[:,delay_index])[max_field_index]
    y1 = y[delaytime_range]

    theory_E_delay = np.array(abs_data[:, delay_index])
    param = model.fit(y1, weights=(theory_E_delay/max(theory_E_delay))[delaytime_range], x=x1, phi_0=0, a=0, b=0.0005, c=0, d=0)

    return x, y, x1, y1, param

def real_time_phase_plot(mat: str, kyvalues: list, delay_index: int, low_lim: float, high_lim: float):
    x, y, x1, y1, param = phase_fitting(mat, kyvalues, delay_index, low_lim, high_lim)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    fig.set_size_inches(10, 8)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # print(np.arange(-200,210,10)[delay_index])
    plt.scatter(x, y, color='red')
    ax.plot(x1, [phase_fun(t, param.best_values['phi_0'], param.best_values['a'], param.best_values['b'], param.best_values['c'], param.best_values['d']) for t in x1])
    # plt.plot(x1, param.eval(x=x1))
    plt.xlabel("Time (fs)")
    plt.ylabel("Phase (radians)")
    plt.xlim([-70,70])
    plt.ylim([-2,2])
    # plt.savefig("phase.png")
    # plt.show()
    return ax

def time_delay_phase_data(mat: str, kyvalues: list, delay_array, low_lim: float, high_lim: float):
    # _, abs_data, _, _ = extract_data(mat, kyvalues)
    # theory_b_parameter = np.zeros(abs_data.shape[1])
    # theory_c_parameter = np.zeros(abs_data.shape[1])
    # delay_array = np.arange(10, -10.2, -0.2)

    def temp_function(delay_index):
        _, _, _, _, param = phase_fitting(mat, kyvalues, delay_index, low_lim, high_lim)
        return param

    # for delay_index in range(len(delay_array)):
    #     _, _, _, _, theory_param = phase_fitting(mat, kyvalues, delay_index, low_lim, high_lim)
    #     theory_b_parameter[delay_index] = theory_param.best_values['b']
    #     theory_c_parameter[delay_index] = theory_param.best_values['c']

    theory_params = Parallel(n_jobs=8)(delayed(temp_function)(i) for i in range(len(delay_array)))

    theory_b_parameter = [obj.best_values['b'] for obj in theory_params]
    theory_c_parameter = [obj.best_values['c'] for obj in theory_params]

    return theory_b_parameter, theory_c_parameter

# def time_delay_phase_plot(mat: str, kyvalues: list, low_lim: float, high_lim: float):
#     fig = plt.figure()
#     ax1 = fig.add_subplot(1, 1, 1)
#     ax2 = fig.add_subplot(1, 1, 1)
#     fig.set_size_inches(10, 8)
#     ax1.spines['left'].set_position('center')
#     ax1.spines['bottom'].set_position('zero')
#     ax1.spines['right'].set_color('none')
#     ax1.spines['top'].set_color('none')
#     ax1.xaxis.set_ticks_position('bottom')
#     ax1.yaxis.set_ticks_position('left')

#     ax2.spines['left'].set_position('center')
#     ax2.spines['bottom'].set_position('zero')
#     ax2.spines['right'].set_color('none')
#     ax2.spines['top'].set_color('none')
#     ax2.xaxis.set_ticks_position('bottom')
#     ax2.yaxis.set_ticks_position('left')

#     # plt.plot(delaytime_array, b_parameter,'o-')
#     theory_b_parameter, theory_c_parameter = time_delay_phase_data(mat, kyvalues, np.arange(10,-10.2,-0.2),low_lim, high_lim)
#     ax1.plot(np.arange(10,-10.2,-0.2), theory_b_parameter, '-*')
#     ax2.plot(np.arange(10,-10.2,-0.2), theory_c_parameter, '-*')
#     # plt.axhline(y=probe_param.best_values['b'], color='r', linestyle='-')
#     plt.legend(['Numerical calculation', 'Probe parameter'], fontsize=20)
#     plt.xlabel("Time delay(fs)", fontsize=25)
#     plt.ylabel("b parameter", fontsize=25)
#     plt.xlim([-10,10])
#     ax1.yaxis.set_label_coords(0, 0.7)
#     ax1.tick_params(axis='x', labelsize=25)
#     ax1.tick_params(axis='y', labelsize=25)
#     ax1.set_yticks([-0.0002, 0.0006])
#     # plt.ylim([-0.0009, 0.0015])
#     # plt.savefig("phase.png")
#     # plt.show()
#     return ax1, ax2
if __name__ == '__main__':
    print(extract_data('MgO', [0]))

# ax = real_time_plot('MgO', [0], 50)
# plt.show()
# print(real_time_plot('MgO', [0], 50))
# fig, ax = plt.subplots()
# ax1 = time_delay_plot('MgO', [0])
# ax2 = time_delay_plot('ZnO', [0])

# fig.axes.append(time_delay_plot('MgO', [0]))
# fig.axes.append(time_delay_plot('ZnO', [0]))
# plt.xlabel("Time Delay (fs)", fontsize=25)
# plt.ylabel("Summed E-field ", fontsize=25)
# ax1.tick_params(axis='x', labelsize=25)
# ax1.tick_params(axis='y', labelsize=25)
# plt.legend(['Numerical Data'],fontsize=18)
# plt.xlim([-10.2,10])
# plt.show()