import numpy as np

# hyperelastic constants
# *Hyperelastic, neo hooke, moduli=instantaneous
c10 = 4.249e2
d1 = 5.043e-4
# Prony parameters (g_i, tau_i)
prony = [
    [8.962e-02, 1.075e+00],
    [1.189e-01, 8.694e+00],
    [1.354e-01, 7.030e+01],
    [1.048e-01, 5.684e+02],
    [9.308e-02, 4.596e+03],
    [6.789e-02, 3.717e+04],
    [4.156e-02, 3.005e+05],
    [4.097e-02, 2.430e+06],
    [3.949e-02, 1.965e+07],
    [4.396e-02, 1.589e+08],
    [4.428e-02, 1.285e+09],
    [4.339e-02, 1.039e+10],
    [4.199e-02, 8.401e+10],
    [3.879e-02, 6.793e+11],
    [2.487e-02, 5.493e+12],
    [1.488e-02, 4.441e+13]]
# TRS parameters (log10(a_T), temperature)
trs = [
    [13.479, -40.06],
    [12.223, -30.00],
    [11.508, -20.02],
    [10.351, -10.02],
    [9.168, -0.01],
    [7.906, 9.99],
    [6.510, 19.94],
    [4.932, 30.10],
    [2.965, 39.98],
    [0.000, 49.99],
    [-2.427, 59.99],
    [-4.238, 69.99],
    [-5.398, 79.94],
    [-6.806, 89.98],
    [-8.595, 100.00],
    [-10.499, 109.99],
    [-12.222, 119.99],
    [-14.000, 130.00],
    [-15.398, 139.94],
    [-16.523, 149.95]]
# lower limit for relaxation times after shift
relaxationTimeLimit = 1.

# files to be created
pronyInputFile = open('parameter-prony.inp', 'wt')
prfInputFile = open('parameter-prf.inp', 'wt')
userInputFile = open('parameter-user.inp', 'wt')

# end of input parameters

prony = np.array(prony)
trs = np.array(trs)
nProps = trs.size + 1

hyper = '*Hyperelastic, neo hooke, moduli=instantaneous\n' + \
        f'{c10:.3e}, {d1:.3e}\n'
pronyInputFile.write(hyper)
prfInputFile.write(hyper)
userInputFile.write(hyper)

pronyInputFile.write('*Viscoelastic, time=prony\n')

g0 = 2 * c10
trsProny = '*TRS, definition=tabular\n'
trsParasUser = []
column = 1
for loga, temp in trs:
    trsProny += f'{loga:.3f}, {temp:.2f}\n'
    trsParasUser += loga, temp
trsUserText = ''
for para in trsParasUser:
    if column % 8 == 0:
        trsUserText += '\n'
    else:
        trsUserText += ', '
    trsUserText += f'{para:.3f}'
    column += 1

i = 0
for gi, taui in prony:
    i += 1
    pronyInputFile.write(f'{gi:.3e},, {taui:.3e}\n')
    a = 1 / taui / 3 / g0 / gi
    nonlinear = f'** relaxation time {taui:.3e}:\n' + \
                '*Viscoelastic, nonlinear, networkid=' + \
                f'{i}, sratio={gi:.3e}, law=strain\n{a:.3e}, 1, 0\n'
    prfInputFile.write(nonlinear + trsProny)
    userInputFile.write(nonlinear)
    minShift = relaxationTimeLimit / taui
    userInputFile.write('*TRS, definition=user, ' + \
                        f'properties={nProps}\n{minShift:.3e}' + \
                        trsUserText + '\n')

pronyInputFile.write(trsProny)

pronyInputFile.close()
prfInputFile.close()
userInputFile.close()
