# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import math
import flow_stress

R = 8.3144598

burgers = 2.5E-10
TAILOR = 3.1
ALPHA_DISL = 0.15

A_D1 = 0.2
A_D2 = 1.2

A_REX = 8.6E-4

A_REX_GG = 0.64
B_REX_GG = 8.3E-4
M_REX_0 = 5.4E-11

Qsd_parameters = {
    "Qsd_main_parametres": {
        'Qsd_intercept': 311691.0,
        'Qsd_exp_coef': -278242.0,
        'Qsd_exp_index': -0.394
    },
    "Qsd_elem_parametres_list": {
        'Mn': [88752.0, 0.310],
        'Si': [22801.0, 1.000],
        'Mo': [84864.0, 0.650],
        'Ni': [-38575.0, 0.300],
        'V': [-7298.0, 1.000],
        'Nb': [132594.0, 0.263],
        'Ti': [82128.0, 0.401],
        'Cr': [0.0, 0.000],
        'W': [0.0, 0.000],
        'Al': [0.0, 0.000],
        'Co': [0.0, 0.000]
    }
}


Qd_params = {
    "intercept": 267000.0,
    "C_coef": -2535.52,
    "Mn_coef":  1010.0,
    "Si_coef":  33620.76,
    "Mo_coef": 35651.28,
    "Ti_coef": 93680.52,
    "Ti_index": 0.5919,
    "V_coef": 31673.46,
    "Nb_index": 0.5649,
    "Nb_coef": 70729
}

chemComposition = {
    "C": 0.15,
    "Mn": 0.55,
    "Si": 0.26,
    "Mo": 0.01,
    "Ti": 0.020,
    "V": 0.030,
    "Nb": 0.020,
    "Zr": 0.03,
    "Al": 0.01,
    "N": 0.006
}
'''
chemComposition = {
    "C": 0.11,
    "Mn": 0.36,
    "Si": 0.23,
    "Cr": 0.39,
    "Mo": 0.31,
    "Cu": 0.46,
    "Ni": 1.9,
    "Nb": 0.01,
    "V": 0.01,
    "N": 0.006
}
'''
A_params = {
    "intercept": 12.197,
    "C": 65.590,
    "Nb": -49.052
}

pic_def_params_n = {
    "intercept": 0.161842770261532,
    "C": -0.00696611705120967,
    "Si": 0.000411126099894321,
    "Mn": -0.00195610427350708,
    "Mo": -0.0414779247496787,
    "Ti": -0.219036103186581,
    "V": -0.336962200415879,
    "Nb": 0.416986308305169,
    "N": 2.41384210803522
}

pic_def_params_A = {
    "intercept": 0.0917675486280388,
    "C": 0.00926549708330233,
    "Si": -0.000911941287290171,
    "Mn": 0.00214928161680221,
    "Mo": 0.069515483004259,
    "Ti": 0.0349197422162702,
    "V": -0.0244100370067512,
    "Nb": -0.0584538744552952,
    "N": 0.441598554601606
}

return_parametrs = {
    "v_D": 2.0E12,
    "mu_pre_coef": 2.6,
    "U_rec": 248100.0,
    "V_rec_factor": 2900.0,
    "V_rec_exp": -5.5E-4
}

parameters = {
    "Qd_params": Qd_params,
    "A_params": A_params,
    "pic_def_params_n": pic_def_params_n,
    "pic_def_params_A": pic_def_params_A,
    "Qsd_parameters": Qsd_parameters,
    "return_parametrs": return_parametrs
}


class Chemical:
    def __init__(self, chem_composition, base="Fe", lattice_type="FCC"):
        self.bale = base
        self.lattice_type = lattice_type
        self.chemComposition = chem_composition
        if base not in self.chemComposition:
            self.chemComposition[base] = self.__calc_base_conc__(chem_composition)
        self.atom_mass = {
            'Fe': 55.847,
            'C': 12.0107,
            'Si': 28.086,
            'Mn': 54.9381,
            'Cr': 51.996,
            'Ni': 58.71,
            'Mo': 95.94,
            'Cu': 63.54,
            'V': 50.942,
            'Nb': 92.906,
            'Ti': 47.90,
            'Zr': 91.22,
            'Co': 58.9331944,
            'O': 15.9994,
            'S': 32.064,
            'P': 30.9737622,
            'H': 1.00784,
            'B': 10.806,
            'N': 14.0067,
            'Al': 26.9815,
            'Ca': 40.08,

        }
        self.mol_concentations, self.MR = self.mol_concentations()
        self.y_concentrations = self.y_concentrations(self.mol_concentations)

    def __calc_base_conc__(self, imput_composition):
        output = 100
        for i in imput_composition:
            output -= imput_composition[i]
        if output < 0:
            output = 0
        return output

    def mol_concentations(self):
        m = 0
        mr = 0
        for i in self.chemComposition:
            if i in self.atom_mass:
                m += self.chemComposition[i] / self.atom_mass[i]
                mr += self.chemComposition[i] * self.atom_mass[i]
        mr /= 100
        mol_conc = {}
        for i in self.chemComposition:
            if i in self.atom_mass:
                mol_conc[i] = (self.chemComposition[i] / self.atom_mass[i]) / m
        return mol_conc, mr

    def y_concentrations(self, mol_concentations):
        a = 0.0
        b = 0.0
        if self.lattice_type == "FCC":
            a = 1.0
            b = 1.0
        elif self.lattice_type == "BCC":
            a = 1.0
            b = 3.0
        elif self.lattice_type == "trigonal":
            a = 1.0
            b = 3.0

        intro_elements = 0
        for i in mol_concentations:
            if i == "C" or i == "N" or i == "H" or i == "B" or i == "S":
                intro_elements += mol_concentations[i]
        y_list = {}
        for i in mol_concentations:
            if i == "C" or i == "N" or i == "H" or i == "B" or i == "S":
                y_el = (a/b)*(mol_concentations[i]/(1 - intro_elements))
                y_list[i] = y_el
                y_list["va"] = 1 - y_el
            else:
                y_el = mol_concentations[i]/(1 - intro_elements)
                y_list[i] = y_el
        return y_list


class QSD:
    def __init__(self, input_parameters, chem_composition, base="Fe", lattice_type="FCC"):
        self.parameters = input_parameters["Qsd_parameters"]
        self.chem = Chemical(chem_composition, base=base, lattice_type=lattice_type)
        self.qsd = self.qsd_austenite_calc()

    def qsd_austenite_calc(self):
        q = self.parameters["Qsd_main_parametres"]["Qsd_intercept"]
        if "C" in self.chem.y_concentrations:
            a = self.parameters["Qsd_main_parametres"]["Qsd_exp_coef"]
            b = self.parameters["Qsd_main_parametres"]["Qsd_exp_index"]
            q += a * (1.0 - math.exp(b * self.chem.y_concentrations["C"]))
        for i in self.chem.y_concentrations:
            if i in self.parameters["Qsd_elem_parametres_list"]:
                q += self.parameters["Qsd_elem_parametres_list"][i][0]*pow(self.chem.y_concentrations[i], self.parameters["Qsd_elem_parametres_list"][i][1])
        return q



class Dislocations:
    def __init__(self, parameters, composition):
        self.parameters = parameters
        self.composition = composition

        self.fs = flow_stress.FlowStress(parameters=parameters, chem_composition=composition)


    def austenite_disl_dens(self, d_s=0, T=1.0, C=0.0):
        #s = self.fs.calcFlowStress(D_0=D_0, e=e, v=v, T=T)
        #s_y = self.fs.calcFlowStress(D_0=D_0, e=0.2, v=v, T=T)
        mu = self.mu(T-273.15)
        a = (self.austenite_a(T=T, C=C), )
        b =  self.burgers(params=a) # (burgers, ) #

        return pow(d_s /(ALPHA_DISL*TAILOR*mu*b[0]), 2.0)

    def mu(self, t):
        return 8.1E10 * (0.91 - ((t - 300)/1810))

    def burgers(self, params=(0.356E-9, 0.356E-9, ), alpha=math.pi/2.0, lattice_type="FCC"):
        try:
            if lattice_type == "FCC":
                 return(0.5*params[0], )
            elif lattice_type == "BCC":
                return ((3.0/4.0)*params[0], params[0], )
            elif lattice_type == "HCP":
                return (params[0], params[1])
            elif lattice_type == "trigonal":
                t = math.cos(alpha / 2.0)
                return (params[0] , 4 * params[0] * t*t, )
            elif lattice_type == "tetragonal":
                return (params[0], params[1])
            elif lattice_type == "simple_cubic":
                return (params[0], )
            else:
                return None
        except:
            return None

    def austenite_a(self, T=1273.15, C=0.0):
        return (0.356 + 0.0000085*T + 0.0033*C)*pow(10, -9)

    def austenite_disl_dens_t(self, x_t=0.0, D_0=20, e=0.8, v=0.0011, T=1000.0, C=0.0, tau=0.0):
        s = self.fs.calcFlowStress(D_0=D_0, e=e, v=v, T=T)
        s_y = self.fs.calcFlowStress(D_0=D_0, e=0.002, v=v, T=T)

        new_ds = self.disl_dens_decr_austenite(s_0=s, d_s=(s - s_y), T=T, C=C, tau=tau)

        #print("ds= ", new_ds)

        ro_avg = self.austenite_disl_dens(d_s=new_ds, T=T, C=C)

        #print("rp_avg= ", ro_avg)

        return A_D1 * ro_avg * (1 + (1.0/A_D1)*math.exp(-A_D2 * x_t))

    def austenite_p_rex(self, x_t=0.0, D_0=20, e=0.8, v=0.0011, T=1000.0, C=0.0, tau=0.0):
        dens = self.austenite_disl_dens_t(x_t=x_t, D_0=D_0, e=e, v=v, T=T, C=C, tau=tau)
        #print("ro_disl ", dens)
        mu = self.mu(T - 273.15)
        a = (self.austenite_a(T=T, C=C),)
        b = self.burgers(params=a)
        return 0.5*mu*b[0]*b[0]*dens

    def n_rex(self, x_t=0.0, D_0=20, e=0.8, v=0.0011, T=1000.0, C=0.0, tau=0.0):
        p_rex = self.austenite_p_rex(x_t=x_t, D_0=D_0, e=e, v=v, T=T, C=C, tau=tau)
        return (A_REX * p_rex * p_rex)/D_0

    def m_rex_gb(self, T):
        q_sd = QSD(self.parameters, self.composition).qsd
        q_rex_gg = A_REX_GG * q_sd
        s_rex_gg = B_REX_GG * q_rex_gg #*2.05
        #print("q_rex_gg ", math.exp((s_rex_gg/R))*math.exp(-q_rex_gg/(R*T)))
        #print("s_rex_gg ", math.exp(-q_rex_gg/(R*T)))
        return M_REX_0*math.exp((s_rex_gg/R))*math.exp(-q_rex_gg/(R*T))

    def disl_dens_decr_austenite(self, s_0=1.0, d_s=0.0, T=1000.0, C=0.1, tau=0.0):
        mu = self.mu(T-273.15)
        v_d = self.parameters["return_parametrs"]["v_D"]
        u_rec = self.parameters["return_parametrs"]["U_rec"]
        a = (self.austenite_a(T=T, C=C),)
        b = self.burgers(params=a)
        v_rec = self.v_rec(T=T, b=b[0])
        e_t = self.parameters["return_parametrs"]["mu_pre_coef"]*mu
        a_p = ALPHA_DISL
        M = TAILOR
        factor_s = (64.0 * pow(d_s, 2) * v_d)/(9.0 * pow(M, 3) * pow(a_p, 2) * e_t)
        #print("a_p ", a_p)
        exp_s = math.exp(-u_rec/(R*T))
        #print("exp_s= ", exp_s)
        #print("u_rec= ", u_rec)
        sin_s = math.sinh((d_s * v_rec)/(R*T))
        #print("R*T= ", R*T)

        return s_0 - (factor_s * exp_s * sin_s * tau)

    def v_rec(self, T=1000.0, b=3.0E-9):
        a = self.parameters["return_parametrs"]["V_rec_factor"]
        e = self.parameters["return_parametrs"]["V_rec_exp"]
        return a * pow(b, 3) * math.exp(e*T)


class RecrystalizationSolver:
    def __init__(self, parameters, chem_composition):
        self.dislocations_handler = Dislocations(parameters, chem_composition)
        self.D_0 = None
        self.e = None
        self.v = None
        self.T = None
        self.is_set = False

    def set_initial(self, D_0=20, e=0.8, v=0.0011, T=1000.0):
        self.D_0 = D_0
        self.e = e
        self.v = v
        self.T = T
        self.is_set = True

    def solve_austenite_plot(self, n_step_max=10000000000000000, d_tau=1.0, max_tau=1000000000000000.0, max_x=0.999):
        x_t = 0
        x_t_prev = 0
        inner_m_p = 0
        tau = 0
        n_step = 1
        #chem = Chemical(self.dislocations_handler.composition)
        o_x = [tau, ]
        o_y = [x_t, ]
        o_y_2 = [None, ]
        o_y_3 = [None, ]
        c = self.dislocations_handler.composition["C"]
        #n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=self.D_0, e=self.e, v=self.v, T=self.T, C=c)
        n_rex = None
        #print(tau, x_t)
        print(x_t)
        dens = self.dislocations_handler.austenite_disl_dens_t(x_t=x_t, D_0=self.D_0, e=self.e, v=self.v,
                                                               T=self.T, C=c, tau=tau)

        d_current = self.D_0

        if self.is_set:
            while n_step < n_step_max and tau < max_tau and x_t < max_x:
                dx_t = 100.0
                p_rex = self.dislocations_handler.austenite_p_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)
                #print(p_rex)
                n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)
                d_current = n_rex**(1.0/3.0)
                #print(n_rex**(1.0/3.0))
                m_rex_gb = self.dislocations_handler.m_rex_gb(self.T)
                #print(m_rex_gb)
                inner_m_p_tmp = inner_m_p
                inner_m_p_tmp += m_rex_gb * p_rex * d_tau
                #print("m_p ", inner_m_p)
                x_ext_t = n_rex * pow(inner_m_p_tmp, 3)
                x_t = 1 - math.exp(-x_ext_t)

                while dx_t > 0.00001:
                    x_t = (x_t + x_t_prev) / 2.0
                    p_rex = self.dislocations_handler.austenite_p_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T,
                                                                      C=c, tau=tau)
                    n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)

                    d_current = n_rex ** (1.0 / 3.0)

                    inner_m_p_tmp = inner_m_p
                    inner_m_p_tmp += m_rex_gb * p_rex * d_tau
                    x_ext_t = n_rex * pow(inner_m_p_tmp, 3)

                    dx_t = abs(x_t - x_t_prev)/x_t
                    x_t_prev = x_t
                    x_t = 1 - math.exp(-x_ext_t)
                print(x_t)
                inner_m_p += m_rex_gb * p_rex * d_tau

                dens = self.dislocations_handler.austenite_disl_dens_t(x_t=x_t, D_0=d_current, e=self.e, v=self.v,
                                                                       T=self.T, C=c, tau=tau)

                tau += d_tau
                n_step += 1
                o_x.append(tau)
                o_y.append(x_t)
                o_y_2.append(n_rex ** (1.0 / 3.0))
                o_y_3.append(dens)
        print(x_t)
        o_x.append(tau)
        o_y.append(x_t)
        o_y_2.append(n_rex ** (1.0 / 3.0))
        o_y_3.append(dens)
        plt.title("plot")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.scatter(o_x, o_y, edgecolors='r', s=10)

        plt.grid(True, linestyle='-', color='0.75')
        plt.show()

        plt.grid(True, linestyle='-', color='0.75')
        plt.scatter(o_x, o_y_2, edgecolors='b', s=10)
        plt.show()

        plt.grid(True, linestyle='-', color='0.75')
        plt.scatter(o_x, o_y_3, edgecolors='g', s=10)
        plt.show()

    def solve_austenite(self, n_step_max=10000000000000000, d_tau=1, max_tau=1000000000000000.0, max_x=0.999):
        x_t = 0
        x_t_prev = 0
        inner_m_p = 0
        tau = 0
        n_step = 1

        x_t_out = {tau: x_t}
        c = self.dislocations_handler.composition["C"]
        #n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=self.D_0, e=self.e, v=self.v, T=self.T, C=c)

        d_current = self.D_0

        if self.is_set:
            while n_step < n_step_max and tau < max_tau and x_t < max_x:
                dx_t = 100.0
                p_rex = self.dislocations_handler.austenite_p_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)
                #print(p_rex)
                n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)
                d_current = n_rex**(1.0/3.0)
                #print(n_rex**(1.0/3.0))
                m_rex_gb = self.dislocations_handler.m_rex_gb(self.T)
                #print(m_rex_gb)
                inner_m_p_tmp = inner_m_p
                inner_m_p_tmp += m_rex_gb * p_rex * d_tau
                #print("m_p ", inner_m_p)
                x_ext_t = n_rex * pow(inner_m_p_tmp, 3)
                x_t = 1 - math.exp(-x_ext_t)

                while dx_t > 0.00001:
                    x_t = (x_t + x_t_prev) / 2.0
                    p_rex = self.dislocations_handler.austenite_p_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T,
                                                                      C=c, tau=tau)
                    n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)

                    d_current = n_rex ** (1.0 / 3.0)

                    inner_m_p_tmp = inner_m_p
                    inner_m_p_tmp += m_rex_gb * p_rex * d_tau
                    x_ext_t = n_rex * pow(inner_m_p_tmp, 3)

                    dx_t = abs(x_t - x_t_prev)/x_t
                    x_t_prev = x_t
                    x_t = 1 - math.exp(-x_ext_t)

                inner_m_p += m_rex_gb * p_rex * d_tau

                dens = self.dislocations_handler.austenite_disl_dens_t(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)
                #print(dens)

                tau += d_tau
                n_step += 1
                x_t_out[tau] = x_t


        x_t_out[tau] = x_t

        return x_t_out

    def x_t_calc(self, x_t, d_current, inner_m_p, c, tau, d_tau, p_z=0.0):

        dx_t = 100.0
        p_rex = self.dislocations_handler.austenite_p_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c,
                                                          tau=tau) - p_z
        # print(p_rex, p_z)
        n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)
        d_current = n_rex ** (1.0 / 3.0)

        m_rex_gb = self.dislocations_handler.m_rex_gb(self.T)
        inner_m_p_tmp = inner_m_p
        inner_m_p_tmp += m_rex_gb * p_rex * d_tau
        x_ext_t = n_rex * pow(inner_m_p_tmp, 3)
        x_t = 1 - math.exp(-x_ext_t)
        x_t_prev = x_t

        while dx_t > 0.00001:
            x_t = (x_t + x_t_prev) / 2.0
            p_rex = self.dislocations_handler.austenite_p_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T,
                                                              C=c, tau=tau)
            n_rex = self.dislocations_handler.n_rex(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T, C=c, tau=tau)

            d_current = n_rex ** (1.0 / 3.0)

            inner_m_p_tmp = inner_m_p
            inner_m_p_tmp += m_rex_gb * p_rex * d_tau
            x_ext_t = n_rex * pow(inner_m_p_tmp, 3)

            if x_t > 0:
                dx_t = abs(x_t - x_t_prev) / x_t
            else:
                dx_t = 0.0
            x_t_prev = x_t
            x_t = 1 - math.exp(-x_ext_t)

        inner_m_p += m_rex_gb * p_rex * d_tau
        dens = self.dislocations_handler.austenite_disl_dens_t(x_t=x_t, D_0=d_current, e=self.e, v=self.v, T=self.T,
                                                               C=c, tau=tau)
        return x_t, inner_m_p, dens, d_current

'''
disl = Dislocations(parameters, chemComposition)

s = disl.fs.calcFlowStress(D_0=20, e=0.8, v=0.0011, T=1000.0)
print(s)
e_p = disl.fs.e_p(D_0=20, v=0.0011, T=1000.0)
print(e_p)
s_y = disl.fs.calcFlowStress(D_0=20, e=0.002, v=0.0011, T=1000.0)
print("s_y=", s_y)
disl_dens = disl.austenite_disl_dens(x_t=0.0, s=s, s_y=s_y, T=1000.0, C=0.1)
print(disl_dens)


n = disl.n_rex(x_t=0, D_0=20, e=0.8, v=0.0011, T=1000.0)
print(n)
n = disl.n_rex(x_t=0.5, D_0=20, e=0.8, v=0.0011, T=1000.0)
print(n)
n = disl.n_rex(x_t=1, D_0=20, e=0.8, v=0.0011, T=1000.0)
print(n)

composition = Chemical(chemComposition)
print(composition.mol_concentations)
print(composition.MR)
print(composition.y_concentrations)

q_sd = QSD(parameters, chemComposition)
print(q_sd.qsd)

m = disl.m_rex_gb(1500.0)
print(m)



solver = RecrystalizationSolver(parameters, chemComposition)
solver.set_initial(D_0=83.0, e=0.35, v=2.0, T=1350)
solver.solve_austenite_plot(d_tau=0.01, max_tau=1000.0)
print(solver.solve_austenite())
'''









