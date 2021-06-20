########################################################################################
# Python script to calculate overpotential and plot step figure of OER and ORR
# Written by Yafei Jiang
# Email: jiangyafei730@163.com
# usage: python OER.py
# directory format: each dir of models should contains all intermediates of ReactionPath
#                   dir of each intermediate should contain opt and freq directories.
########################################################################################
import sys
from math import exp, log
import os
import re
import numpy as np
import matplotlib.pyplot as plt


def calc_Gibbs(filename, T, state):
    '''read OUTCAR of frequency calculation and Temperature, 
       return Ezpe, Hcorr, Gcorr'''
    ###############################################################################
    # calculate free energy contribution (ZPE included) of vibrational modes
    # Modified script from vib2G.py Written by Hai Xiao <haixiao@wag.caltech.edu>
    ###############################################################################

    T = float(T)
    outcar = open(filename)
    state = int(state)  # 0 for stationary point; 1 for transition state

    # Get vibrational eigenvalues from OUTCAR
    vibmatch = "f  ="
    ivibmatch = "f/i="
    Evib = []
    Nivib = 0
    Eivib = []
    for line in outcar:
        if vibmatch in line:
            # Read in vibrational mode in unit of meV (hv, NOT 1/2*hv)
            Evib.append(float(line.split()[9]))
        if ivibmatch in line:
            Eivib.append(float(line.split()[8]))
            Nivib = Nivib + 1
    outcar.close()

    # Calculate characteristic vibrational temperatures
    # conversion factor from http://physics.nist.gov/cgi-bin/cuu/Value?evk
    # using the source: 2014 CODATA recommended values
    # T = hv/k
    meV2K = 11.6045221
    Tvib = []
    for i in Evib:
        Tvib.append(i * meV2K)

    # Set the threshold for low frequency modes, which might lead to unphysically large entropy contribution
    # Currently use the peak frequency 60 cm-1 of ice Ih spectra, corresponding to
    # the acoustic translational mode of the six-member rings
    # Planck constant in eV*s from http://physics.nist.gov/cgi-bin/cuu/Value?hev
    # using the source: 2014 CODATA recommended values
    h = float("4.135667662E-15")
    # speed of light in vacuum in unit of cm/s (exact)
    c = float("2.99792458E10")
    # conversion factor for wavenumber 1/lambda (cm-1): h*c
    cm2eV = h * c
    ###print "%10.6E"%(cm2eV)
    threshold = 50.0 * cm2eV * 1000.0 * meV2K
    ###print "%12.6f"%(threshold)
    for i in range(len(Tvib)):
        if Tvib[i] < threshold:
            print(
                "Warning: very small frequency found, which might lead to unphysically large entropy contribution")
            print("very small frequency %3i: %12.6f K" % (i + 1, Tvib[i]))
            Tvib[i] = threshold
            print("         treated as real: %12.6f K" % (Tvib[i]))

    # Important check of small imaginary frequencies
    if state == 0:
        if Nivib > 0:
            print("Warning: imaginary frequencies found for supposed local minimum !!!")
            print("Warning: all imaginary frequencies treated as threshold ones for now")
            for i in range(Nivib):
                print("imaginary frequency %i: %12.6f K" % (i + 1, Eivib[i] * meV2K))
                Tvib.append(threshold)
                print("       treated as real: %12.6f K" % (threshold))
    elif state == 1:
        if Nivib == 0:
            print("Error: no imaginary frequency found for supposed transition state !!!")
            sys.exit(0)
        elif Nivib > 1:
            print("Warning: more than 1 imaginary frequency found for supposed transition state !!!")
            print("Warning: all imaginary frequencies, except the lowest, treated as threshold ones for now")
            for i in range(Nivib - 1):
                print("imaginary frequency %i: %12.6f K" % (i + 1, Eivib[i] * meV2K))
                Tvib.append(threshold)
                print("       treated as real: %12.6f K" % (threshold))
            print("imaginary frequency %i: %12.6f K" % (Nivib, Eivib[Nivib - 1] * meV2K))
            print("treated as the only 1 imaginary frequency for supposed transition state")
    else:
        print("Error: only accept state 0 or 1 for now")
        sys.exit(0)

    # Free energy contribution, using formulas from
    # http://www.gaussian.com/g_whitepap/thermo.htm
    # R = N_A * k
    # use k here, in order to get final energies in the unit of eV
    # Boltzmann constant from http://physics.nist.gov/cgi-bin/cuu/Value?tkev
    # using the source: 2014 CODATA recommended values
    k = float("8.6173303E-5")
    # ZPE
    H_zpe = 0.0
    for i in Tvib:
        H_zpe = H_zpe + k * i * 0.5

    # enthalpy and entropy contributions
    H_vib = 0.0
    S_vib = 0.0
    if T == 0.0:
        H_vib = 0.0
        S_vib = 0.0
    else:
        for i in Tvib:
            q = i / T
            H_vib = H_vib + k * i / (exp(q) - 1.0)
            S_vib = S_vib + k * (q / (exp(q) - 1.0) - log(1.0 - exp(-q)))

    G_vib = H_vib - T * S_vib
    total_vib = H_zpe + G_vib
    print("ZPE contribution is %20.16f eV" % (H_zpe))
    print("Enthalpy contribution H at %.2f K is %20.16f eV" % (T, H_vib))
    print("Entropy S at %.2f K is %20.16f eV/K" % (T, S_vib))
    print("H - TS at %.2f K is %20.16f eV" % (T, G_vib))
    print("Total contribution at %.2f K with ZPE is %20.16f eV" % (T, total_vib))
    return H_zpe, H_vib, total_vib


def get_Gibbs(ReactionPath):
    '''get Gibbs free energy of all the species at current directory'''
    models = input("Please input file path of models: \n such as: model1 model2 model3 \n Press Enter to skip input models.\n").split()
    if not models:  # if no input, script will scan the current directory and child directory
        # models = []
        for root, dirs, files in os.walk("./"):
            for dir in dirs:
                if dir == ReactionPath[-1]:
                    models.append(root[2:])  # Get names of models
                else:
                    continue
    else:
        pass
    E0_models = []
    Ezpe_models = []
    Hcorr_models = []
    Gcorr_models = []
    G_models = []
    # ReactionPath = ["v","OH","O","OOH"]
    for name in models:
        E0_model = []
        Ezpe_model = []
        Hcorr_model = []
        Gcorr_model = []
        G_model = []
        for species in ReactionPath:
            species_path = name + "/" + species
            if os.path.exists(species_path + "/opt/OUTCAR"):
                file = open(species_path + "/opt/OUTCAR")
                lines = []
                key = "entropy"
                for line in file.readlines():
                    if key in line:
                        lines.append(line)
                s = re.findall(r".?\d+\.?\d*", lines[-1])
                E0_species = float(s[-1])
                file.close()
            else:
                print("Warning: %s is not exist!" %species_path)
                sys.exit()
    
            if os.path.exists(species_path + "/freq/OUTCAR"):
                print(name + "_" + species)
                Ezpe_species, Hcorr_species, Gcorr_species = calc_Gibbs(species_path + "/freq/OUTCAR", Temperature, 0)
            else:
                Ezpe_species, Hcorr_species, Gcorr_species = 0.0, 0.0, 0.0

            E0_model.append(E0_species)
            Ezpe_model.append(Ezpe_species)
            Hcorr_model.append(Hcorr_species)
            Gcorr_model.append(Gcorr_species)
            G_model.append(E0_species + Gcorr_species)
            print("=" * 60)
        
        G_models.append(G_model)
        E0_models.append(E0_model)
        Ezpe_models.append(Ezpe_model)
        Hcorr_models.append(Hcorr_model)
        Gcorr_models.append(Gcorr_model)

    with open("Energy.dat", 'w') as f:
        f.write("{:^10s}{:^14s}{:^14s}{:^14s}{:^14s}{:^14s}\n".format("Species", "E0", "Ezpe", "Hcorr", "Gcorr", "Gtotal"))
        for name in models:
            i = models.index(name)
            f.write("%s \n" %name)
            for species in ReactionPath:
                j = ReactionPath.index(species)
                f.write("{:^10s}{:^14.6f}{:^14.6f}{:^14.6f}{:^14.6f}{:^14.6f}\n".format(species, E0_models[i][j], Ezpe_models[i][j], Hcorr_models[i][j], Gcorr_models[i][j], G_models[i][j]))
    
    return models,np.array(G_models)  # the Gibbs energy of species for all the models


def get_OER(models,ReactionPath,G,pH):
    '''calculate relative Gibbs free energy of species for OER and obtain overpotential'''
    # ReactionPath = ["v","OH","O","OOH"]
    deltaG = -4.92  # eV  2H2 + O2 --> 2H2O
    GOH = GH2O - GH2 / 2
    # GO2 = GH2O * 2 - GH2 * 2 - deltaG  # eV
    G_U0_models = []  # Relatvie Gibbs free energies at U=0V
    G_U123_models = []  # Relatvie Gibbs free energies at U=1.23V
    # G_Ueta_models = []  # Relatvie Gibbs free energies at U=(1.23 + eta) V
    eta_models = []  # Overpotential
    for name in models:  # for each model
        if pH <= 7:
            # acid condition
            deltaG_pH = pH * kb * Temperature * log(10)
            i = models.index(name)
            G_U0_1 = G[i][1] + GH2 / 2 - G[i][0] - GH2O + deltaG_pH  # v + H2O --> *OH + H+ + e-
            G_U0_2 = G[i][2] + GH2 - G[i][0] - GH2O + 2*deltaG_pH  # v + H2O --> *O + 2H+ + 2e-
            G_U0_3 = G[i][3] + GH2 * 3 / 2 - G[i][0] - GH2O * 2 + 3*deltaG_pH  # v + 2H2O --> *OOH + 3H+ + 3e-
            G_U0_4 = -deltaG + 4 * deltaG_pH  # v + 2H2O --> v + O2 + 4H+ + 4e-
            G_U123_1 = G_U0_1 + 1 * deltaG / 4
            G_U123_2 = G_U0_2 + 2 * deltaG / 4
            G_U123_3 = G_U0_3 + 3 * deltaG / 4
            G_U123_4 = G_U0_4 + 4 * deltaG / 4
            eta = max(G_U123_1,G_U123_2-G_U123_1,G_U123_3-G_U123_2,G_U123_4-G_U123_3)
            print("The overpotential of model %s is %.2f V" %(models[i],eta))
            G_U0_models.append(np.array([0.0,G_U0_1,G_U0_2,G_U0_3,G_U0_4]))
            G_U123_models.append(np.array([0.0,G_U123_1,G_U123_2,G_U123_3,G_U123_4]))
            eta_models.append(eta)
        else: 
            # alkaline condition
            deltaG_pH = (14-pH) * kb * Temperature * log(10)
            i = models.index(name)
            G_U0_1 = G[i][1] - G[i][0] - GOH - deltaG_pH  # v + OH- --> *OH + e-
            G_U0_2 = G[i][2] + GH2O - G[i][0] - 2* GOH - 2*deltaG_pH  # v + 2OH- --> *O + H2O + 2e-
            G_U0_3 = G[i][3] + GH2O - G[i][0] - 3* GOH - 3*deltaG_pH  # v + 3OH- --> *OOH + H2O + 3e-
            G_U0_4 = -deltaG - 4 * deltaG_pH  # v + 4OH- --> v + O2 + 2H2O + 4e-
            G_U123_1 = G_U0_1 + 1 * deltaG / 4
            G_U123_2 = G_U0_2 + 2 * deltaG / 4
            G_U123_3 = G_U0_3 + 3 * deltaG / 4
            G_U123_4 = G_U0_4 + 4 * deltaG / 4
            eta = max(G_U123_1,G_U123_2-G_U123_1,G_U123_3-G_U123_2,G_U123_4-G_U123_3)
            print("The overpotential of model %s is %.2f V" %(models[i],eta))
            G_U0_models.append(np.array([0.0,G_U0_1,G_U0_2,G_U0_3,G_U0_4]))
            G_U123_models.append(np.array([0.0,G_U123_1,G_U123_2,G_U123_3,G_U123_4]))
            eta_models.append(eta)
    
    G_U0 = np.array(G_U0_models)  # Relatvie Gibbs free energies at U=0V
    G_U123 = np.array(G_U123_models)  # Relatvie Gibbs free energies at U=1.23V
    etas = np.array(eta_models)  # overpotential
    
    write_file(models, ReactionPath[:], G_U0_models,G_U123_models,eta_models)
            
    return G_U0,G_U123


def get_ORR(models,ReactionPath,G,pH):
    '''calculate relative Gibbs free energy of species for OER and obtain overpotential'''
    # ReactionPath = ["v","OOH","O","OH"]
    deltaG = -4.92  # eV  2H2 + O2 --> 2H2O
    GOH = GH2O - GH2 / 2
    GO2 = GH2O * 2 - GH2 * 2 - deltaG  # eV
    G_U0_models = []  # Relatvie Gibbs free energies at U=0V
    G_U123_models = []  # Relatvie Gibbs free energies at U=1.23V
    # G_Ueta_models = []  # Relatvie Gibbs free energies at U=(1.23 + eta) V
    eta_models = []  # Overpotential
    for name in models:  # for each model
        if pH <= 7:
            # acid condition
            deltaG_pH = pH * kb * Temperature * log(10)
            i = models.index(name)
            G_U0_1 = G[i][1] - G[i][0] - GO2 - GH2 / 2 - deltaG_pH  # v + O2 + H+ + e- --> *OOH
            G_U0_2 = G[i][2] + GH2O - G[i][0] - GO2 - GH2 - 2*deltaG_pH  # v + O2 + 2H+ + 2e- --> *O + H2O
            G_U0_3 = G[i][3] + GH2O - G[i][0] - GO2 - GH2 * 3 / 2 - 3*deltaG_pH  # v + O2 + 3H+ + 3e- --> *OH + H2O
            G_U0_4 = deltaG - 4 * deltaG_pH  # v + O2 + 4H+ + 4e- --> v + 2H2O
            G_U123_1 = G_U0_1 - deltaG / 4
            G_U123_2 = G_U0_2 - 2 * deltaG / 4
            G_U123_3 = G_U0_3 - 3 * deltaG / 4
            G_U123_4 = G_U0_4 - 4 * deltaG / 4
            eta = max(G_U123_1,G_U123_2-G_U123_1,G_U123_3-G_U123_2,G_U123_4-G_U123_3)
            print("The overpotential of model %s is %.2f V" %(models[i],eta))
            G_U0_models.append(np.array([0.0,G_U0_1,G_U0_2,G_U0_3,G_U0_4]))
            G_U123_models.append(np.array([0.0,G_U123_1,G_U123_2,G_U123_3,G_U123_4]))
            eta_models.append(eta)
        else: 
            # alkaline condition
            deltaG_pH = (14-pH) * kb * Temperature * log(10)
            i = models.index(name)
            G_U0_1 = G[i][1] + GOH - G[i][0] - GH2O - GO2 + deltaG_pH  # v + O2 + H2O + e- --> *OOH + OH-
            G_U0_2 = G[i][2] + 2 * GOH - G[i][0] - GH2O - GO2 + 2*deltaG_pH  # v + O2 + H2O + 2e- --> *O + 2OH-
            G_U0_3 = G[i][3] + 3 * GOH - G[i][0] - 2 * GH2O - GO2 + 3*deltaG_pH  # v + O2 + 2H2O + 3e- --> *OH + 3OH-
            G_U0_4 = deltaG + 4 * deltaG_pH  # v + O2 + 2H2O + 4e- --> v + 4OH-
            G_U123_1 = G_U0_1 - deltaG / 4
            G_U123_2 = G_U0_2 - 2 * deltaG / 4
            G_U123_3 = G_U0_3 - 3 * deltaG / 4
            G_U123_4 = G_U0_4 - 4 * deltaG / 4
            eta = max(G_U123_1,G_U123_2-G_U123_1,G_U123_3-G_U123_2,G_U123_4-G_U123_3)
            print("The overpotential of model %s is %.2f V" %(models[i],eta))
            G_U0_models.append(np.array([0.0,G_U0_1,G_U0_2,G_U0_3,G_U0_4]))
            G_U123_models.append(np.array([0.0,G_U123_1,G_U123_2,G_U123_3,G_U123_4]))
            eta_models.append(eta)
    
    G_U0 = np.array(G_U0_models)  # Relatvie Gibbs free energies at U=0V
    G_U123 = np.array(G_U123_models)  # Relatvie Gibbs free energies at U=1.23V
    etas = np.array(eta_models)  # overpotential
    
    write_file(models, ReactionPath[:], G_U0_models,G_U123_models,eta_models)
            
    return G_U0,G_U123


def get_NRR(models,ReactionPath,G,pH):
    '''calculate relative Gibbs free energy of species for NRR and obtain overpotential'''
    # ReactionPath = ["v","NN","NNH","NNH2","N","NH","NH2","NH3"]
    deltaG = -0.34  # eV  N2 + 3H2 --> 2NH3
    GN2 = GNH3 * 2 - GH2 * 3 - deltaG  # eV
    G_U0_models = []  # Relatvie Gibbs free energies at U=0V
    G_UV_models = []  # Relatvie Gibbs free energies at U=0.34V
    # G_Ueta_models = []  # Relatvie Gibbs free energies at U=(0.34 + eta) V
    eta_models = []  # Overpotential
    for name in models:  # for each model
        # acid condition
        deltaG_pH = pH * kb * Temperature * log(10)
        i = models.index(name)
        G_U0_1 = G[i][1] - GN2 - G[i][0]  # v + N2 --> *NN
        G_U0_2 = G[i][2] - GN2 - G[i][0] - GH2 / 2 - deltaG_pH  # v + N2 + (H+ + e-) --> *NNH
        G_U0_3 = G[i][3] - GN2 - G[i][0] - GH2 - 2 * deltaG_pH  # v + N2 + 2(H+ + e-) --> *NNH2
        G_U0_4 = G[i][4] - GN2 - G[i][0] - GH2 * 3/2 - 3 * deltaG_pH  # v + N2 + 3(H+ + e-) --> *N + NH3
        G_U0_5 = G[i][5] - GN2 - G[i][0] - GH2 * 2 - 4 * deltaG_pH  # v + N2 + 4(H+ + e-) --> *NH + NH3
        G_U0_6 = G[i][6] - GN2 - G[i][0] - GH2 * 5/2 - 5 * deltaG_pH  # v + N2 + 5(H+ + e-) --> *NH2 + NH3
        G_U0_7 = G[i][7] - GN2 - G[i][0] - GH2 * 3 - 6 * deltaG_pH  # v + N2 + 6(H+ + e-) --> *NH3 + NH3
        G_U0_8 = deltaG - 6*deltaG_pH  # v + N2 + 6(H+ + e-) --> v + 2NH3
        G_UV_1 = G_U0_1
        G_UV_2 = G_U0_2 - 1 * deltaG / 6
        G_UV_3 = G_U0_3 - 2 * deltaG / 6
        G_UV_4 = G_U0_4 - 3 * deltaG / 6
        G_UV_5 = G_U0_5 - 4 * deltaG / 6
        G_UV_6 = G_U0_6 - 5 * deltaG / 6
        G_UV_7 = G_U0_7 - 6 * deltaG / 6
        G_UV_8 = G_U0_8 - 6 * deltaG / 6
        eta = max(G_UV_2-G_UV_1,G_UV_3-G_UV_2,G_UV_4-G_UV_3,G_UV_5-G_UV_4,G_UV_6-G_UV_5,G_UV_7-G_UV_6)
        print("The overpotential of model %s is %.2f V" %(models[i],eta))
        G_U0_models.append(np.array([0.0,G_U0_1,G_U0_2,G_U0_3,G_U0_4,G_U0_5,G_U0_6,G_U0_7,G_U0_8]))
        G_UV_models.append(np.array([0.0,G_UV_1,G_UV_2,G_UV_3,G_UV_4,G_UV_5,G_UV_6,G_UV_7,G_UV_8]))
        eta_models.append(eta)
            
    G_U0 = np.array(G_U0_models)  # Relatvie Gibbs free energies at U=0V
    G_UV = np.array(G_UV_models)  # Relatvie Gibbs free energies at U=1.23V
    etas = np.array(eta_models)  # overpotential
    
    write_file(models, ReactionPath[:], G_U0_models,G_UV_models,eta_models)
            
    return G_U0,G_UV


def write_file(models, ReactionPath, G_U0_models,G_UV_models,eta_models):
    ReactionPath.append(ReactionPath[0])
    with open("Energy.dat", 'a') as f:
        f.write("\n" + "="*40 + "\n")
        f.write("{:^10s}{:^14s}{:^14s}\n".format("Species", "U=0V", "U=1.23V"))
        for name in models:
            i = models.index(name)
            f.write("%s \n" %name)
            for j in range(len(ReactionPath)):
                f.write("{:^10s}{:^14.6f}{:^14.6f}\n".format(ReactionPath[j], G_U0_models[i][j], G_UV_models[i][j]))
        f.write("\n")
        for name in models:
            f.write("The overpotential of model %s is %.2f V\n" %(name,eta_models[models.index(name)]))
        
        
def plot_line_dot(x, y, color, path_label, TextLabel, FontSize=22):
    """绘制虚实折线图"""
    y_max, y_min = np.max(y), np.min(y)  # 获取y值的最大值和最小值
    y_bias = (y_max - y_min) / 50  # 获取文本标签y方向偏移量

    for i in range(len(y)):  # 遍历所有列的能量值
        y_new = []
        x_new = []
        # 1.生成新的XY坐标点，个数加倍
        for j in range(len(y[i])):
            y_new.append(y[i][j])
            y_new.append(y[i][j])
            x_new.append(2*j+1)
            x_new.append(2*j+2)
        # 2.绘制实线折线图
        k = 0
        while k < len(y_new):
            x_line = [x_new[k], x_new[k+1]]
            y_line = [y_new[k], y_new[k+1]]
            plt.plot(x_line, y_line, linestyle='-', linewidth=4, color=color[i])
            k += 2

        # 3.绘制虚线折线图
        plt.plot(x_new, y_new, linestyle='--', linewidth=2, color=color[i], label=path_label[i])
        # 4.标记能量值，偏移量视具体情况而定
        if TextLabel:
            for j in range(len(x)):
                plt.text(x[j] * 2 - 0.9, y[i][j] + y_bias, "{:.2f}".format(y[i][j]), fontsize=FontSize, color=color[i])


def plot_Reaction(models,ReactionPath,y):
    x = [i+1 for i in range(len(ReactionPath)+1)]
    ReactionPath[0] = ""
    ReactionPath.append("")
    ReactionPath = ["*" + s for s in ReactionPath]

    plt.figure(figsize=(12, 8), dpi=300)  # 设置图片大小及分辨率
    
    TextLabel = False  # 是否添加坐标点对应的数值文本标签
    FontSize = 20  # 设置能量值文本大小
    Axis_FontSize = 18  # 设置XY轴标签及标题大小
    color = ["red","purple","orange","green","cyan","blue","brown","pink","yellow","black"]  # add new colors if models > 10
    path_label = models
    plot_line_dot(x, y, color, path_label, TextLabel, FontSize)

    plt.xlim(x[0] * 2 - 1.5, x[-1] * 2 + 1)  # x轴刻度范围
    plt.xticks([i * 2 - 0.5 for i in x], ReactionPath, fontsize=Axis_FontSize)  # x轴标签
    plt.yticks(fontsize=Axis_FontSize)  # y轴标签
    y_max, y_min = np.max(y), np.min(y)  # 获取y值的最大值和最小值
    y_scale = (y_max - y_min) / 10  # y轴延伸长度
    plt.ylim(y_min - y_scale, y_max + y_scale)  # y轴刻度范围    
    plt.ylabel("Relative Gibbs free energy (eV)", fontsize=Axis_FontSize)  # 纵轴标题
    # plt.xlabel(X_title, fontsize=Axis_FontSize)  # 横轴标题
    # plt.title(pic_title, fontsize=Axis_FontSize)  # 图标题
    
    plt.legend(fontsize=Axis_FontSize-2, loc="lower right")  # 添加图例,位置在右下角
    # plt.legend(fontsize=Axis_FontSize-2, loc="center left",bbox_to_anchor=(-0.08,-0.2),ncol=2)  # 图例放在图外
    plt.tight_layout()
    
    #plt.show()  # 展示图片
    plt.savefig("EnergyProfile.png")  # 保存图片到当前目录


def main():
    Reaction = input("Please input the type of eletrochemical reaction: \n [1] OER \n [2] ORR \n [3] NRR \n ")
    pH = float(input("Please input the pH value: \n"))
    if Reaction == "1":
        # ReactionPath = input("Please input the intermediates of the Reaction: \n such as: v OH O OOH \n").split()
        ReactionPath = ["v","OH","O","OOH"]
        models,G = get_Gibbs(ReactionPath)
        G_U0,G_UV = get_OER(models,ReactionPath,G,pH)
        plot_Reaction(models,ReactionPath[:],G_UV)
    elif Reaction == "2":
        # ReactionPath = input("Please input the intermediates of the Reaction: \n such as: v OOH O OH \n").split()
        ReactionPath = ["v","OOH","O","OH"]
        models,G = get_Gibbs(ReactionPath)
        G_U0,G_UV = get_ORR(models,ReactionPath,G,pH)
        plot_Reaction(models,ReactionPath[:],G_UV)
    elif Reaction == "3":
        ReactionPath = input("Please input the intermediates of the Reaction: \n such as: slab NN NNH NNH2 N NH NH2 NH3 \n").split()
        models,G = get_Gibbs(ReactionPath)
        G_U0,G_UV = get_NRR(models,ReactionPath,G,pH)
        plot_Reaction(models,ReactionPath[:],G_UV)
    else:
        print("The type of reaction you input is not supported. \n")
        sys.exit()
        
        
if __name__ == '__main__':
    global kb, Temperature, GH2, GH2O, GNH3
    # use kb here, in order to get final energies in the unit of eV
    # Boltzmann constant from http://physics.nist.gov/cgi-bin/cuu/Value?tkev
    # using the source: 2014 CODATA recommended values
    kb = float("8.6173303E-5")  # eV/K
    Temperature = 298.15  # K
    # The Gibbs free energies of H2,H2O and NH3 depend on your computational methods.
    GH2 = -6.805  # eV
    GH2O = -14.229  # eV at 298.15K and 0.035bar
    GNH3 = -19.10  # eV  
    
    main()