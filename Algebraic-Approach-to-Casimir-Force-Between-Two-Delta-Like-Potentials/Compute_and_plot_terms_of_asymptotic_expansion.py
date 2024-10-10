#!/usr/bin/python3


##############################
# The script for computing and drawing terms of asymptotic expanisons
# of Casimir energy, equation (67) in paper
# https://doi.org/10.1007/s00023-020-00994-2

# Computation were made in Python 3.10.12???
##############################


##############################
import numpy as np
# import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
##############################



# print("numpy version: ", np.__version__)
# print("scipy version: ", scipy.__version__)



def compute_and_plot_1(gamma0, gammaN, numberOfGammaPoints,
                       yMinPlot, yMaxPlot, nameOfPlot, nameOfErrorPlot):
    """Function computes and plot values of four terms of asymptotic
    expansion presented in the equations (67) as the functions of
    parameter gamma.

    Terms 2*alpha/pi^3, chi/lambda and chi*b_1 are ignored in the
    computations.
    """


    # Values that parameter gamma will take.
    gamma_array = np.linspace(gamma0, gammaN, numberOfGammaPoints)

    # Array for storing the results of computations.
    # Results of evaluation first term
    firstTermArray = np.zeros(numberOfGammaPoints)
    errorFirstTermArray = np.zeros(numberOfGammaPoints)

    secondTermArray = np.zeros(numberOfGammaPoints)
    errorSecondTermArray = np.zeros(numberOfGammaPoints)

    thirdTermArray = np.zeros(numberOfGammaPoints)
    errorThirdTermArray = np.zeros(numberOfGammaPoints)

    fourthTermArray = np.zeros(numberOfGammaPoints)
    errorFourthTermArray = np.zeros(numberOfGammaPoints)

    # Line below is only a reminder to the reader.
    # resultAndError = (0, 0)



    # Loop for computing values of functions for various values of gamma
    for i, gamma in enumerate(gamma_array):
        # First function to integrate
        I_1 = lambda l: \
            np.exp(-2 * l) / ((gamma + l) *
                              ((gamma + l)**2 - np.exp(-2 * l)))

        # Integrating of first function.
        resultAndError = integrate.quad(I_1, 0, np.inf)

        firstTermArray[i] = resultAndError[0]
        errorFirstTermArray[i] = resultAndError[1]



        # Second function to integrate
        I_2 = lambda l: \
            (l**2) * (3 * (gamma + l)**2 * np.exp(-2 * l) -
                      np.exp(-4 * l)) / ((gamma + l)**2 *
                                         ((gamma + l)**2 -
                                          np.exp(-2 * l))**2)

        # Integrating of second function.
        resultAndError = integrate.quad(I_2, 0, np.inf)

        secondTermArray[i] = resultAndError[0] / gamma
        errorSecondTermArray[i] = resultAndError[1] / gamma



        # Third function to integrate
        I_3 = lambda l: \
            l * np.exp(-2 * l) / ((gamma + l) * ((gamma + l)**2 -
                                                 np.exp(-2 * l)))

        # Integrating of third function.
        resultAndError = integrate.quad(I_3, 0, np.inf)

        thirdTermArray[i] = 2 * resultAndError[0] / gamma
        errorThirdTermArray[i] = 2* resultAndError[1] / gamma



        # Fourth function to integrate
        I_4 = lambda l: \
            (1 - l) * np.exp(-2 * l) / ((gamma + l)**2 -
                                        np.exp(-2 * l))

        resultAndError = integrate.quad(I_4, 0, np.inf)

        fourthTermArray[i] = resultAndError[0] / gamma
        errorFourthTermArray[i] = resultAndError[1] / gamma



    fig, axis1 = plt.subplots()

    axis1.spines[["right", "top"]].set_visible(False)

    axis1.plot(gamma_array, firstTermArray, 'r', label=r"$I_{ 1 }$")
    axis1.plot(gamma_array, secondTermArray, 'g', label=r"$I_{ 2 }$")
    axis1.plot(gamma_array, thirdTermArray, 'b', label=r"$I_{ 3 }$")
    axis1.plot(gamma_array, fourthTermArray, 'y', label=r"$I_{ 4 }$")

    axis1.set_ylim(yMinPlot, yMaxPlot)

    # EN:
    descriptionOfXAxis = r"$\gamma$, range [{}, {}]".format(gamma0,
                                                            gammaN)
    # PL:
    # descriptionOfXAxis = r"$\gamma$, zakres [{}, {}]".format(gamma0,
    #                                                          gammaN)

    axis1.set_xlabel(descriptionOfXAxis)

    # EN:
    axis1.set_ylabel("Numerical values")
    # PL:
    # axis1.set_ylabel("Wartości numeryczne")

    axis1.legend()

    fig.savefig(nameOfPlot)

    axis1.cla()



    axis1.plot(gamma_array, errorFirstTermArray, 'r',
             label=r"$\Delta I_{ 1 }$")
    axis1.plot(gamma_array, errorSecondTermArray, 'g',
             label=r"$\Delta I_{ 2 }$")
    axis1.plot(gamma_array, errorThirdTermArray, 'b',
             label=r"$\Delta I_{ 3 }$")
    axis1.plot(gamma_array, errorFourthTermArray, 'y',
             label=r"$\Delta I_{ 4 }$")

    axis1.set_xlabel(descriptionOfXAxis)

    # EN:
    # axis1.set_ylabel("Values of numerical errors")
    # PL:
    plt.ylabel("Wartości numeryczne błędów")

    axis1.legend()

    fig.savefig(nameOfErrorPlot)

    plt.close()





def compute_and_plot_2(gamma0, gammaN, numberOfGammaPoints,
                       yMinPlot, yMaxPlot, nameOfPlot):
    """Function computes and plot sum of four terms of asymptotic
    expansion presented in the equations (67) as the functions of
    parameter gamma.

    Terms 2*alpha/pi^3 is ignored. We choose arbitrary values chi/lambda = 10
    and chi*b_1 = 1. This last two values are consisten with
    inequalities derived in PhD thesis (only in polish) "Algebraiczne
    podejście do efektu Casimira w układzie z dwoma potencjałami typu
    quasi delta" by Kamil Ziemian.
    """


    # Values that parameter gamma will take.
    gamma_array = np.linspace(gamma0, gammaN, numberOfGammaPoints)

    # Array for storing the results of computations.
    # Results of evaluation first term
    firstTermArray = np.zeros(numberOfGammaPoints)
    errorFirstTermArray = np.zeros(numberOfGammaPoints)

    secondTermArray = np.zeros(numberOfGammaPoints)
    errorSecondTermArray = np.zeros(numberOfGammaPoints)

    thirdTermArray = np.zeros(numberOfGammaPoints)
    errorThirdTermArray = np.zeros(numberOfGammaPoints)

    fourthTermArray = np.zeros(numberOfGammaPoints)
    errorFourthTermArray = np.zeros(numberOfGammaPoints)

    # Line below is only a reminder to the reader.
    # resultAndError = (0, 0)



    # Loop for computing values of functions for various values of gamma
    for i, gamma in enumerate(gamma_array):
        # First function to integrate
        I_1 = lambda l: \
            np.exp(-2 * l) / ((gamma + l) *
                              ((gamma + l)**2 - np.exp(-2 * l)))

        # Integrating of first function.
        resultAndError = integrate.quad(I_1, 0, np.inf)

        firstTermArray[i] = resultAndError[0]
        errorFirstTermArray[i] = resultAndError[1]



        # Second function to integrate
        I_2 = lambda l: \
            (l**2) * (3 * (gamma + l)**2 * np.exp(-2 * l) -
                      np.exp(-4 * l)) / ((gamma + l)**2 *
                                         ((gamma + l)**2 -
                                          np.exp(-2 * l))**2)

        # Integrating of second function.
        resultAndError = integrate.quad(I_2, 0, np.inf)

        secondTermArray[i] = resultAndError[0] / gamma
        errorSecondTermArray[i] = resultAndError[1] / gamma



        # Third function to integrate
        I_3 = lambda l: \
            l * np.exp(-2 * l) / ((gamma + l) * ((gamma + l)**2 -
                                                 np.exp(-2 * l)))

        # Integrating of third function.
        resultAndError = integrate.quad(I_3, 0, np.inf)

        thirdTermArray[i] = 2 * resultAndError[0] / gamma
        errorThirdTermArray[i] = 2* resultAndError[1] / gamma



        # Fourth function to integrate
        I_4 = lambda l: \
            (1 - l) * np.exp(-2 * l) / ((gamma + l)**2 -
                                        np.exp(-2 * l))

        resultAndError = integrate.quad(I_4, 0, np.inf)

        fourthTermArray[i] = resultAndError[0] / gamma
        errorFourthTermArray[i] = resultAndError[1] / gamma



    fig, axis1 = plt.subplots()

    axis1.spines[["right", "top"]].set_visible(False)

    Casimir_energy_array = (10.0 * firstTermArray + secondTermArray
        - thirdTermArray + fourthTermArray)

    axis1.plot(gamma_array, Casimir_energy_array, 'b',
               label=r"$E_{ \text{Cas} }$")
    # axis1.plot(gamma_array, secondTermArray, 'g', label=r"$I_{ 2 }$")
    # axis1.plot(gamma_array, thirdTermArray, 'b', label=r"$I_{ 3 }$")
    # axis1.plot(gamma_array, fourthTermArray, 'y', label=r"$I_{ 4 }$")

    axis1.set_ylim(yMinPlot, yMaxPlot)

    # EN:
    descriptionOfXAxis = (r"$\gamma$," + " range [{}, {}], ".format(gamma0,
                                                                    gammaN)
                          + r" $\frac{ \chi }{ \lambda } = 10.0$, $\chi b_{ 1 } = 1.0$")
    # PL:
    # descriptionOfXAxis = r"$\gamma$, zakres [{}, {}]".format(gamma0,
    #                                                          gammaN)

    axis1.set_xlabel(descriptionOfXAxis)

    # EN:
    axis1.set_ylabel("Numerical values")
    # PL:
    # axis1.set_ylabel("Wartości numeryczne")

    axis1.legend()

    fig.savefig(nameOfPlot)

    axis1.cla()



    # axis1.plot(gamma_array, errorFirstTermArray, 'r',
    #          label=r"$\Delta I_{ 1 }$")
    # axis1.plot(gamma_array, errorSecondTermArray, 'g',
    #          label=r"$\Delta I_{ 2 }$")
    # axis1.plot(gamma_array, errorThirdTermArray, 'b',
    #          label=r"$\Delta I_{ 3 }$")
    # axis1.plot(gamma_array, errorFourthTermArray, 'y',
    #          label=r"$\Delta I_{ 4 }$")

    axis1.set_xlabel(descriptionOfXAxis)

    # EN:
    # axis1.set_ylabel("Values of numerical errors")
    # PL:
    # plt.ylabel("Wartości numeryczne błędów")

    # axis1.legend()

    # fig.savefig(nameOfErrorPlot)

    plt.close()





def compute_first_function_value(gamma):
    """Function computes value of the first term in asymptotic expansion
    (67) at choosen value of parameter gamma.

    Terms 2*alpha/pi^3 and chi/lambda are ignored in the computations.
    """

    I_1 = lambda l: \
        np.exp(-2 * l) / ((gamma + l) *
                          ((gamma + l)**2 - np.exp(-2 * l)))

    return integrate.quad(I_1, 0, np.inf)




def compute_and_plot_sum_1(gamma0, gammaN, numberOfGammaPoints,
                           yMinPlot, yMaxPlot, nameOfPlot):
    """Function computes and plot values of sum of third and fourth
    term of asymptotic expansion presented in the equations (67) as
    the functions of parameter gamma.

    Term 2*alpha/pi^3 is ignored in the computations.
    """


    # Values that parameter gamma will take.
    gamma_array = np.linspace(gamma0, gammaN, numberOfGammaPoints)

    thirdTermArray = np.zeros(numberOfGammaPoints)

    fourthTermArray = np.zeros(numberOfGammaPoints)

    # Line below is only a reminder to the reader.
    # sumArray = np.zeros(numberOfGammaPoints)

    # Line below is only a reminder to the reader.
    # resultAndError = (0, 0)



    # Loop for computing values of functions for various values of gamma
    for i, gamma in enumerate(gamma_array):
        # Third function to integrate
        I_3 = lambda l: \
            l * np.exp(-2 * l) / ((gamma + l) * ((gamma + l)**2 -
                                                 np.exp(-2 * l)))

        # Integrating of third function.
        resultAndError = integrate.quad(I_3, 0, np.inf)

        thirdTermArray[i] = 2 * resultAndError[0] / gamma



        # Fourth function to integrate
        I_4 = lambda l: \
            (1 - l) * np.exp(-2 * l) / ((gamma + l)**2 -
                                        np.exp(-2 * l))

        resultAndError = integrate.quad(I_4, 0, np.inf)

        fourthTermArray[i] = resultAndError[0] / gamma



    sumArray = -thirdTermArray + fourthTermArray

    fig, axis1 = plt.subplots()

    axis1.spines[["right", "top"]].set_visible(False)


    axis1.plot(gamma_array, sumArray, 'b', label=r"$-I_{ 3 } + I_{ 4 }$")

    axis1.set_ylim(yMinPlot, yMaxPlot)

    # EN:
    # descriptionOfXAxis = r"$\gamma$, range [{}, {}]".format(gamma0,
    #                                                         gammaN)
    # PL:
    descriptionOfXAxis = r"$\gamma$, zakres [{}, {}]".format(gamma0,
                                                             gammaN)

    axis1.set_xlabel(descriptionOfXAxis)

    # EN:
    plt.ylabel("Numerical values")
    # PL:
    # axis1.set_ylabel("Wartości numeryczne")

    axis1.legend()

    fig.savefig(nameOfPlot)

    axis1.cla()





def compute_and_plot_ratio_1(gamma0, gammaN, numberOfGammaPoints,
                             yPlotMin, yPlotMax, nameOfPlot):
    """Function computes and plot ratio of first and second term of
    asymptotic expansion presented in the equations (67).

    Terms 2*alpha/pi^3, chi/lambda and chi*b_1 are ignored in the
    computations.
    """


    # Values that parameter gamma will take.
    gamma_array = np.linspace(gamma0, gammaN, numberOfGammaPoints)

    # Array for storing the results of computations contains results
    # Results of evaluation first term
    firstTermArray = np.zeros(numberOfGammaPoints)
    errorFirstTermArray = np.zeros(numberOfGammaPoints)

    secondTermArray = np.zeros(numberOfGammaPoints)
    errorSecondTermArray = np.zeros(numberOfGammaPoints)



    # Loop for computing values of functions for various values of gamma
    for i, gamma in enumerate(gamma_array):
        # First function to integrate
        I_1 = lambda l: \
            np.exp(-2 * l) / ((gamma + l) *
                              ((gamma + l)**2 - np.exp(-2 * l)))

        # Integrating of first function.
        resultAndError = integrate.quad(I_1, 0, np.inf)

        firstTermArray[i] = resultAndError[0]


        # Second function to integrate
        I_2 = lambda l: \
            (l**2) * (3 * (gamma + l)**2 * np.exp(-2 * l) -
                      np.exp(-4 * l)) / ((gamma + l)**2 *
                                         ((gamma + l)**2 -
                                          np.exp(-2 * l))**2)

        # Integrating of second function.
        resultAndError = integrate.quad(I_2, 0, np.inf)

        secondTermArray[i] = resultAndError[0] / gamma



    ratioExpressionArray = firstTermArray / secondTermArray


    plt.plot(gamma_array, ratioExpressionArray, 'b',
             label=r"$\frac{ I_{ 1 } }{ I_{ 2 } }$")

    # EN:
    # descriptionOfXAxis = r"$\gamma$, range [{}, {}]".format(gamma0,
    #                                                         gammaN)
    # PL:
    descriptionOfXAxis = r"$\gamma$, zakres [{}, {}]".format(gamma0,
                                                             gammaN)

    plt.xlabel(descriptionOfXAxis)

    # EN:
    plt.ylabel("Numerical values")
    # PL:
    # plt.ylabel("Wartości numeryczne")

    plt.ylim(yPlotMin, yPlotMax)


    plt.legend()

    plt.savefig(nameOfPlot)

    plt.close()





# print(compute_first_function_value(1.01))
# value: 0.93


gamma0 = 1.1
gammaN = 7.0
numberOfGammaPoints = 1401

compute_and_plot_1(gamma0, gammaN, numberOfGammaPoints, 0, 0.4,
                   "Terms_of_asymptotic_expansion_01.png",
                   "Terms_of_asymptotic_expansion_errors_01.png")

compute_and_plot_2(gamma0, gammaN, numberOfGammaPoints, 0, 5,
                   "Casimir_energy_asymptotic_expansion_01.png")

compute_and_plot_sum_1(gamma0, gammaN, numberOfGammaPoints, 0, 0.25,
                       "Sum_of_third_and_fourth_term_01.png")



gamma0 = 1.01
gammaN = 1.1
numberOfGammaPoints = 200

compute_and_plot_1(gamma0, gammaN, numberOfGammaPoints, 0, 0.95,
                   "Terms_of_asymptotic_expansion_02.png",
                   "Terms_of_asymptotic_expansion_errors_02.png")

compute_and_plot_sum_1(gamma0, gammaN, numberOfGammaPoints, 0, 0.75,
                       "Sum_of_third_and_fourth_term_02.png")



gamma0 = 1.01
gammaN = 7.0
numberOfGammaPoints = 1201

compute_and_plot_ratio_1(gamma0, gammaN, numberOfGammaPoints, 0, 50,
                         "Ratio_of_first_and_second_term.png")

# Ration min: 7.12.
