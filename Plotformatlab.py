
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager

if __name__ == "__main__":
    font = font_manager.FontProperties(family='Times New Roman',
                                       style='normal', size=22)


    origdata = np.loadtxt('origdata.csv', delimiter=',')
    predata = np.loadtxt('predata.csv', delimiter=',')
    traindata = np.loadtxt('traindata.csv', delimiter=',')

    sig = np.loadtxt('sigdata.csv', delimiter=',')

    upf = predata[:, 1] + sig
    lowf = predata[:, 1] - sig


    plt.plot(origdata[:, 0], origdata[:, 1], label='Objective')
    plt.plot(predata[:, 0], predata[:, 1], label='Prediction')
    plt.scatter(traindata[:, 0], traindata[:, 1], 80, c ='none', linewidths=2, edgecolors='k', label='Training data')
    plt.fill_between(origdata[:, 0], upf.ravel(), lowf.ravel(), color='orange', alpha=0.5)
    plt.legend(loc='upper center', prop=font, framealpha= 0.3)
    plt.xlabel('X', fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    h = plt.ylabel('F', fontsize=20)
    h.set_rotation(0)



    plt.show()