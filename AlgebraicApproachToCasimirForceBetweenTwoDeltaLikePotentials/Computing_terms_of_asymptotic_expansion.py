#!/usr/bin/python3


##############################
# Computing and drawing terms of asymptotic expanisons of energy.
# Equation (67) in https://doi.org/10.1007/s00023-020-00994-2

# Computation were made for Python 3????
##############################


##############################
import numpy as np
# import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
##############################



# print("numpy version: ", np.__version__)
# print("scipy version: ", scipy.__version__)



def compute_and_plot_I(gamma0, gammaN, numberOfGammaPoints, nameOfPlot):
    # num_chart = number of created chart
    # alpha = 5.0  # "Coupling" constant of delta potential
    # prob_const = 1.0 / alpha
    # Problematic constant: integral of M_{ p } / lambda

    # prob_const_al = prob_const / alpha

    # N = 100  # Number of points on plot
    # gamma_0 = 1.1  # Lower limit of gamma values. Always must be true gamma > 1.
    # gamma_N = 10.0  # Upper limit of gamma values

    # Values that parameter gamma will take.
    gamma_array = np.linspace(gamma0, gammaN, numberOfGammaPoints)

    # Array for storing the results of computations contains results
    # Results of evaluation first term
    firstTermArray = np.zeros(numberOfGammaPoints)
    secondTermArray = np.zeros(numberOfGammaPoints)
    thirdTermArray = np.zeros(numberOfGammaPoints)
    fourthTermArray = np.zeros(numberOfGammaPoints)

    # Method quad returns tuple so this is fiting.
    resultAndError = (0, 0)



    # Loop for computing integrals for various values of gamma
    for i, gamma in enumerate(gamma_array):
        # First function to integrate
        fun_I = lambda l: \
            np.exp(-2 * l) / ((gamma + l) *
                              ((gamma + l)**2 - np.exp(-2 * l)))

        # Integrating of first function.
        resultAndError = integrate.quad(fun_I, 0, np.inf)

        firstTermArray[i] =  resultAndError[0]



        # Second function to integrate
        fun_II = lambda l: \
            (l**2) * (3 * (gamma + l)**2 * np.exp(-2 * l) -
                      np.exp(-4 * l)) / ((gamma + l)**2 *
                                         ((gamma + l)**2 -
                                          np.exp(-2 * l))**2)

        # Integrating of second function.
        resultAndError = integrate.quad(fun_II, 0, np.inf)

        secondTermArray[i] = resultAndError[0] / gamma



        # Third function to integrate
        fun_III = lambda l: \
            l * np.exp(-2 * l) / ((gamma + l) * ((gamma + l)**2 -
                                                 np.exp(-2 * l)))

        # Integrating of third function.
        resultAndError = integrate.quad(fun_III, 0, np.inf)

        thirdTermArray[i] = -2 * resultAndError[0] / gamma



        # Fourth function to integrate
        fun_IV = lambda l: \
            (1 - l) * np.exp(-2 * l) / ((gamma + l)**2 -
                                        np.exp(-2 * l))

        resultAndError = integrate.quad(fun_IV, 0, np.inf)

        fourthTermArray[i] = resultAndError[0] / gamma



    plt.plot(gamma_array, firstTermArray, 'r', label=r"$f_{ 1 }$")
    plt.plot(gamma_array, secondTermArray, 'g', label=r"$f_{ 2 }$")
    plt.plot(gamma_array, thirdTermArray, 'b', label=r"$f_{ 3 }$")
    plt.plot(gamma_array, fourthTermArray, 'y', label=r"$f_{ 4 }$")

    plt.xlabel(r"$\gamma$")
    plt.ylabel("Skala porównawcza")
    # plt.ylabel(r"Energy [cm$^{ -1 }$]")

    plt.legend()

    # plt.subplot(3, 1, 2)
    # plt.plot(gamma_array, result_list_II)
    # plt.title("Second term")
    # plt.xlabel(r"$\gamma$ [dimensionless]")
    # plt.ylabel(r"Energy [cm$^{ -1 }$]")
    # plt.ylim(0, 0.01)

    # plt.subplot(3, 1, 3)
    # plt.plot(gamma_array, result_list_III)
    # plt.title("Third term")
    # plt.xlabel(r"$\gamma$ [dimensionless]")
    # plt.ylabel(r"Energy [cm$^{ -1 }$]")


    # ax_1.plot(gamma_array, result_list)
    # ax_1.set_title("Total energy")
    # ax_1.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_1.set_ylabel(r"Energy [cm$^{ -1 }$]")

    # ax_2.plot(gamma_array, result_list_I)
    # ax_2.set_title("Energy of first term")
    # ax_2.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_2.set_ylabel(r"Energy [cm$^{ -1 }$]")

    # ax_3.plot(gamma_array, result_list_II)
    # ax_3.set_title("Energy of second term")
    # ax_3.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_3.set_ylabel(r"Energy [cm$^{ -1 }$]")
    # # ax_3.axis('tight')

    # ax_4.plot(gamma_array, result_list_III)
    # ax_4.set_title("Energy of third term")
    # ax_4.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_4.set_ylabel(r"Energy [cm$^{ -1 }$]")
    # # ax_4.axis('tight')

    # fig.tight_layout()
    # plt.subplots_adjust(top=0.87, left=0.18, hspace=0.90)
    # fig.axis('tight')

    plt.savefig("Asymptotic_expansions_terms.png")

    plt.close()


def compute_and_plot_II(gamma0, gammaN, numberOfGammaPoints):
    # num_chart = number of created chart
    # alpha = 5.0  # "Coupling" constant of delta potential
    # prob_const = 1.0 / alpha
    # Problematic constant: integral of M_{ p } / lambda

    # prob_const_al = prob_const / alpha

    # N = 100  # Number of points on plot
    # gamma_0 = 1.1  # Lower limit of gamma values. Always must be true gamma > 1.
    # gamma_N = 10.0  # Upper limit of gamma values

    # Values that parameter gamma will take.
    gamma_array = np.linspace(gamma0, gammaN, numberOfGammaPoints)

    # Array for storing the results of computations contains results
    # Results of evaluation first term
    firstTermArray = np.zeros(numberOfGammaPoints)
    funArray = np.zeros(numberOfGammaPoints)
    # secondTermArray = np.zeros(numberOfGammaPoints)
    # thirdTermArray = np.zeros(numberOfGammaPoints)
    # fourthTermArray = np.zeros(numberOfGammaPoints)

    # Method quad returns tuple so this is fiting.
    resultAndError = (0, 0)



    # Loop for computing integrals for various values of gamma
    for i, gamma in enumerate(gamma_array):
        # First function to integrate
        fun_I = lambda l: \
            np.exp(-2 * l) / ((gamma + l) *
                              ((gamma + l)**2 - np.exp(-2 * l)))

        # Integrating of first function.
        resultAndError = integrate.quad(fun_I, 0, np.inf)

        firstTermArray[i] =  resultAndError[0]


        # Function 1 / gamma^3
        funArray[i] = 1 / (gamma**3)



        # # Second function to integrate
        # fun_II = lambda l: \
        #          (l**2) * (3 * (gamma + l)**2 * np.exp(-2 * l) -
        #                    np.exp(-4 * l)) / ((gamma + l)**2 *
        #                                       ((gamma + l)**2 -
        #                                        np.exp(-2 * l))**2)

        # # Integrating of second function.
        # resultAndError = integrate.quad(fun_II, 0, np.inf)

        # secondTermArray[i] = resultAndError[0] / gamma



        # # Third function to integrate
        # fun_III = lambda l: \
        #          l * np.exp(-2 * l) / ((gamma + l) * ((gamma + l)**2 -
        #                                               np.exp(-2 * l)))

        # # Integrating of third function.
        # resultAndError = integrate.quad(fun_III, 0, np.inf)

        # thirdTermArray[i] = -2 * resultAndError[0] / gamma



        # # Fourth function to integrate
        # fun_IV = lambda l: \
        #           (1 - l) * np.exp(-2 * l) / ((gamma + l)**2 -
        #                                       np.exp(-2 * l))

        # resultAndError = integrate.quad(fun_IV, 0, np.inf)

        # fourthTermArray[i] = resultAndError[0] / gamma



    plt.plot(gamma_array, firstTermArray, 'r', label=r"$f_{ 1 }$")
    plt.plot(gamma_array, funArray, 'g',
             label=r"$\frac{ 1 }{ \gamma^{ 4 } }$")
    # plt.plot(gamma_array, thirdTermArray, 'b', label=r"$f_{ 3 }$")
    # plt.plot(gamma_array, fourthTermArray, 'y', label=r"$f_{ 4 }$")

    plt.xlabel(r"$\gamma$")
    plt.ylabel("Skala porównawcza")
    # plt.ylabel(r"Energy [cm$^{ -1 }$]")

    plt.legend()

    # plt.subplot(3, 1, 2)
    # plt.plot(gamma_array, result_list_II)
    # plt.title("Second term")
    # plt.xlabel(r"$\gamma$ [dimensionless]")
    # plt.ylabel(r"Energy [cm$^{ -1 }$]")
    # plt.ylim(0, 0.01)

    # plt.subplot(3, 1, 3)
    # plt.plot(gamma_array, result_list_III)
    # plt.title("Third term")
    # plt.xlabel(r"$\gamma$ [dimensionless]")
    # plt.ylabel(r"Energy [cm$^{ -1 }$]")


    # ax_1.plot(gamma_array, result_list)
    # ax_1.set_title("Total energy")
    # ax_1.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_1.set_ylabel(r"Energy [cm$^{ -1 }$]")

    # ax_2.plot(gamma_array, result_list_I)
    # ax_2.set_title("Energy of first term")
    # ax_2.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_2.set_ylabel(r"Energy [cm$^{ -1 }$]")

    # ax_3.plot(gamma_array, result_list_II)
    # ax_3.set_title("Energy of second term")
    # ax_3.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_3.set_ylabel(r"Energy [cm$^{ -1 }$]")
    # # ax_3.axis('tight')

    # ax_4.plot(gamma_array, result_list_III)
    # ax_4.set_title("Energy of third term")
    # ax_4.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_4.set_ylabel(r"Energy [cm$^{ -1 }$]")
    # # ax_4.axis('tight')

    # fig.tight_layout()
    # plt.subplots_adjust(top=0.87, left=0.18, hspace=0.90)
    # fig.axis('tight')

    plt.savefig("Asymptotic_expansions_terms.png")

    plt.close()


def compute_and_draw_II(alpha, prob_const, gamma_0, gamma_N, N_gamma,
                     num_chart):
    # num_chart = number of created chart
    # alpha = 5.0  # "Coupling" constant of delta potential
    # prob_const = 1.0 / alpha
    # Problematic constant: integral of M_{ p } / lambda

    prob_const_al = prob_const / alpha

    # N = 100  # Number of points on plot
    # gamma_0 = 1.1  # Lower limit of gamma values. Always must be true gamma > 1.
    # gamma_N = 10.0  # Upper limit of gamma values


    gamma_array = np.linspace(gamma_0, gamma_N, N_gamma)

    # Array to contains results
    # result_list = np.zeros(N)
    result_list_I = np.zeros(N_gamma)  # Results of intergrating first term
    result_list_II = np.zeros(N_gamma)
    result_list_III = np.zeros(N_gamma)
    # error_list = np.zeros(N)
    error_list_I = np.zeros(N_gamma)  # Error of intergrating first term
    error_list_II = np.zeros(N_gamma)
    error_list_III = np.zeros(N_gamma)
    # Czy całkowity błąd to suma błędów?

    # one_ov_a_array = np.zeros(n)

    # Some numbers for computations
    al_two_pi_sq = alpha / (2 * (pi**2))
    res_err = (0, 0)  # quad returns tuple so this is fiting
    mult_const = 1.0  # Constant to multiplu results
    # result = 0.0
    # error = 0.0


    # Loop making true computations
    for i, gamma in enumerate(gamma_array):
        # First function to integrate
        fun_I = lambda t: \
            exp(-2 * t) / ((gamma + t) * ((gamma + t)**2 - exp(-2 * t)))

        # result part I
        res_err = integrate.quad(fun_I, 0, inf)

        mult_const = prob_const / (2 * (pi**2))
        result_list_I[i] =  mult_const * res_err[0]
        error_list_I[i] = mult_const * res_err[1]


        # Second function to integrate
        fun_II = lambda t: \
            t * exp(-2 * t) / ((gamma + t) * ((gamma + t)**2
                                               - exp(-2 * t)))

        # result part II
        res_err = integrate.quad(fun_II, 0, inf)

        mult_const = al_two_pi_sq * 1.0 / (gamma * 2 * (pi**2))
        result, error =  mult_const * res_err[0], mult_const * res_err[1]

        result_list_II[i] = mult_const * res_err[0]
        error_list_II[i] = mult_const * res_err[1]


        # Third function to integrate
        fun_III = lambda t: \
            (t - 1) * exp(-2 * t) / ((gamma + t)**2 - exp(-2 * t))

        res_err = integrate.quad(fun_III, 0, inf)

        mult_const = al_two_pi_sq * 1.0 / (4 * pi * gamma)

        result_list_III[i] = mult_const * res_err[0]
        error_list_III[i] = mult_const * res_err[1]


    # Adding results
    result_list = result_list_I + result_list_II + result_list_III

    # fig, ((ax_1, ax_2), (ax_3, ax_4)) = plt.subplots(2, 2, sharex='col',
    #                                                  sharey='row')

    # # Ploting 4 charts
    title_all = r"Results for $\alpha = {:.1f}$, term $\int M_p\, dp / ( \alpha \lambda ) = {:.1f}$".format(alpha, prob_const / alpha)

    plt.suptitle(title_all)

    plt.subplot(2, 1, 1)
    plt.plot(gamma_array, result_list)
    plt.title("Total energy")
    plt.xlabel(r"$\gamma$ [dimensionless]")
    plt.ylabel(r"Energy [cm$^{ -1 }$]")

    plt.subplot(2, 1, 2)
    plt.plot(gamma_array, result_list_I, 'g', label="I term")
    plt.plot(gamma_array, result_list_II, 'y', label="II term")
    plt.plot(gamma_array, result_list_III, 'b', label="III term")
    plt.plot(gamma_array, result_list, 'r', label="Sum of terms")
    plt.title("Contribution of terms")
    plt.xlabel(r"$\gamma$ [dimensionless]")
    plt.ylabel(r"Energy [cm$^{ -1 }$]")
    plt.legend()
    # plt.ylim(0, 0.01)

    # plt.subplot(3, 1, 3)
    # plt.plot(gamma_array, result_list_III)
    # plt.title("Third term")
    # plt.xlabel(r"$\gamma$ [dimensionless]")
    # plt.ylabel(r"Energy [cm$^{ -1 }$]")


    # ax_1.plot(gamma_array, result_list)
    # ax_1.set_title("Total energy")
    # ax_1.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_1.set_ylabel(r"Energy [cm$^{ -1 }$]")

    # ax_2.plot(gamma_array, result_list_I)
    # ax_2.set_title("Energy of first term")
    # ax_2.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_2.set_ylabel(r"Energy [cm$^{ -1 }$]")

    # ax_3.plot(gamma_array, result_list_II)
    # ax_3.set_title("Energy of second term")
    # ax_3.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_3.set_ylabel(r"Energy [cm$^{ -1 }$]")
    # # ax_3.axis('tight')

    # ax_4.plot(gamma_array, result_list_III)
    # ax_4.set_title("Energy of third term")
    # ax_4.set_xlabel(r"$\gamma$ [dimensionless]")
    # ax_4.set_ylabel(r"Energy [cm$^{ -1 }$]")
    # # ax_4.axis('tight')

    # fig.tight_layout()
    plt.subplots_adjust(top=0.87, left=0.15, hspace=0.60)
    # fig.axis('tight')

    plt.savefig("Casimir_Energ_good_{}.png".format(
        str(num_chart).zfill(2)))
    # plt.close(fig)


def compute_for_probl_const(alpha, prob_const_array, gamma_0, gamma_N,
                            N_gamma, num_chart_series):

    for j, prob_const in enumerate(prob_const_array):
        print(num_chart_series, j, num_chart_series * 5 + j + 1)
        compute_and_draw(alpha, prob_const, gamma_0, gamma_N, N_gamma,
                         num_chart_series * 5 + j + 1)


gamma0 = 1.1
gammaN = 10.0
numberOfGammaPoints = 2001

compute_and_plot_I(gamma0, gammaN, numberOfGammaPoints)
# compute_and_plot_II(gamma0, gammaN, numberOfGammaPoints)
