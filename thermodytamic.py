# -*- coding: utf-8 -*-

import math
import copy
import operator
import rec_press
import functools

R = 8.3144598
epsilon = 0.000001
max_iterations = 20
ponamarenko_const = 0.008314


class CarbonitrideThermodynamicsInSteel():
    def __init__(self, chem_composition={}, td_params={}, t=1000.0):
        self.chem_composition = chem_composition
        if "C" not in self.chem_composition:
            self.chem_composition["C"] = 0
        if "N" not in self.chem_composition:
            self.chem_composition["N"] = 0
        if "Ti" not in self.chem_composition:
            self.chem_composition["Ti"] = 0
        if "Nb" not in self.chem_composition:
            self.chem_composition["Nb"] = 0
        if "V" not in self.chem_composition:
            self.chem_composition["V"] = 0
        if "Zr" not in self.chem_composition:
            self.chem_composition["Zr"] = 0
        self.t = t
        self.chem = rec_press.Chemical(self.chem_composition)
        self.td_params = td_params
        self.excange_energies = self.calculate_excange_energies()

        if t < 1200:
            tmp = - 0.000041117591 * (self.t ** 2) + 0.048518757316 * self.t - 12.499747647620
            self.convergention_ctiterial = math.exp(tmp)
        elif t < 1450:
            tmp = -0.00000017 * (self.t ** 3) + 0.00053087 * (self.t ** 2) - 0.57055987 * self.t + 200.94623558
            self.convergention_ctiterial = math.exp(tmp)
        else:
            self.convergention_ctiterial = 1e-11

    @functools.lru_cache(maxsize=128, typed=False)
    def calc_solute_energy(self, element="C", t=1000.0):
        if "solute_energies" in self.td_params:
            if element in self.td_params["solute_energies"]:
                intercept = self.td_params["solute_energies"][element]["intercept"]
                coef = self.td_params["solute_energies"][element]["coef"]
                return intercept + coef * t
        return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def calc_formation_energy(self, compound="Ti_C", t=1000.0):
        energy = 0
        if "formation_energies" in self.td_params:
            if compound in self.td_params["formation_energies"]:
                for key in self.td_params["formation_energies"][compound]:
                    energy += self.td_params["formation_energies"][compound][key] * (pow(t, float(key)))
        return energy

    @functools.lru_cache(maxsize=128, typed=False)
    def calc_smelting_energy(self, element="Ti", t=1000.0):
        if "melting_temperatures" in self.td_params and "h_melting" in self.td_params:
            if element in self.td_params["melting_temperatures"] and element in self.td_params["h_melting"]:
                return (self.td_params["h_melting"][element]
                       - t*self.td_params["h_melting"][element]/self.td_params["melting_temperatures"][element])
        return 0.0

    @functools.lru_cache(maxsize=128, typed=False)
    def calc_full_energy(self, compound="Ti_C", t=1000.0):
        solute_energy = math.fsum([self.calc_solute_energy(i, t) for i in compound.split("_")])
        formation_energy = self.calc_formation_energy(compound=compound, t=t)
        return formation_energy - solute_energy #- smelting_energy

    def c_n_check(self):
        if self.chem_composition["C"] == 0 and self.chem_composition["N"] == 0:
            return False
        return True

    def check_carbonitride_formers(self):
        if (self.chem_composition["Ti"] == 0 and self.chem_composition["Nb"] == 0 and
                    self.chem_composition["V"] == 0 and self.chem_composition["Zr"] == 0):
            return False
        return True

    @functools.lru_cache(maxsize=128, typed=False)
    def eq_constant(self, composition=""):
        if composition in self.td_params["eq_constant_params"]:
            q = self.td_params["eq_constant_params"][composition]["Q"]
            a = self.td_params["eq_constant_params"][composition]["A"]
            s = self.td_params["eq_constant_params"][composition]["S"]
            return pow(10.0, ((q / self.t) + a * (self.t) - s))
        return 0.0

    @functools.lru_cache(maxsize=128, typed=False)
    def calculate_intaraction_parameter(self, el_i="", el_j="", el_k=None, order=1):
        order_marker = "e"
        if order == 1:
            order_marker = "e"
        elif order == 2:
            order_marker = "r"
        elif order == -2:
            order_marker = "cross"
        else:
            return 0.0

        if el_i in self.td_params["vagner_interactions"]:
            if el_j in self.td_params["vagner_interactions"][el_i][order_marker]:
                val = 0.0
                if "q" in self.td_params["vagner_interactions"][el_i][order_marker][el_j]:
                    val += self.td_params["vagner_interactions"][el_i][order_marker][el_j]["q"] / self.t
                if "s" in self.td_params["vagner_interactions"][el_i][order_marker][el_j]:
                    val += self.td_params["vagner_interactions"][el_i][order_marker][el_j]["s"]
                if el_k is not None and el_k in self.td_params["vagner_interactions"][el_i][order_marker][el_j]:
                    val = self.td_params["vagner_interactions"][el_i][order_marker][el_j][el_k]
                return val
        return 0.0

    def calculate_excange_energies(self):
        output = {}
        for i in self.chem_composition:
            output[i] = {}
            for j in self.chem_composition:
                if i != j:
                    if ((i in self.td_params["energetic_ponomarenko_params"]) and
                        (j in self.td_params["energetic_ponomarenko_params"])):
                        par_1 = self.td_params["energetic_ponomarenko_params"][i]
                        par_2 = self.td_params["energetic_ponomarenko_params"][j]
                        output[i][j] = (1.0 / 2.0) * (pow(((pow(par_1, 0.5)) - (pow(par_2, 0.5))), 2))
        return output

    def activity_coef(self, element="", t=1000.0, max_order=2):
        coef = 1873.0 / t
        sum = 0
        for n in range(max_order):
            for j in self.chem_composition:
                sum += math.pow(self.chem_composition[j], n+1) * self.calculate_intaraction_parameter(el_i=element, el_j=j, order=n+1)

        for j in self.chem_composition:
            for k in self.chem_composition:
                sum += (self.chem_composition[j] * self.chem_composition[k] *
                       self.calculate_intaraction_parameter(el_i=element, el_j=j, el_k=k, order=-2))

        return math.pow(10.0, coef * sum)

    def max_element_consumption(self, el_1="", el_2="C", n_1=1.0, n_2=1.0):
        if el_1 in self.chem_composition and el_2 in self.chem_composition:
            m_1 = self.chem_composition[el_1]
            m_2 = self.chem_composition[el_2]
            if el_1 in self.chem.atom_mass and el_2 in self.chem.atom_mass:
                mr_1 = self.chem.atom_mass[el_1]
                mr_2 = self.chem.atom_mass[el_2]
                m_el_max = (n_1 / n_2) * (m_2 / mr_2) * mr_1
                if m_el_max > m_1:
                    m_el_max = m_1
                return  m_el_max
            else:
                return 0.0
        else:
            return 0.0

    def m_i_calc(self, element="C", compounds={}):
        m_i = 0
        for i in compounds:
            if element in i.split("_"):
                mr_marker = i.split("_")[1]
                if mr_marker == element:
                    mr_marker = i.split("_")[0]
                m_i += compounds[i] / self.chem.atom_mass[mr_marker]

        m_i *= self.chem.atom_mass[element]

        return m_i

    def m_compound_calk(self, element="C", compounds={}, m_i=0):
        m_compound = 0
        for i in compounds:
            if element in i.split("_"):
                m_compound += compounds[i]
        m_compound += m_i
        return m_compound

    def find_diff(self, el_1="Ti", el_2="C", compounds={}, compounds_activity={}):
        compoumd_marker = el_1 + "_" + el_2
        if compoumd_marker in compounds:
            m_elem_icompound = compounds[compoumd_marker]

            m_el_icarb = 0.0
            if el_1 + "_C" in compounds:
                m_el_icarb = compounds[el_1 + "_C"]

            m_el_initr = 0.0
            if el_1 + "_N" in compounds:
                m_el_initr = compounds[el_1 + "_N"]

            mr_1 = 0.0
            if el_1 in self.chem.atom_mass:
                mr_1 = self.chem.atom_mass[el_1]
            else:
                self.chem.atom_mass[el_1] = 0.0

            mr_2 = 0
            if el_2 in self.chem.atom_mass:
                mr_2 = self.chem.atom_mass[el_2]
            else:
                self.chem.atom_mass[el_2] = 0.0

            mci = self.m_i_calc(element="C", compounds=compounds)
            mni = self.m_i_calc(element="N", compounds=compounds)
            m_i = 0

            m_1_n = 0.0
            if el_1 in self.chem_composition:
                m_1_n = self.chem_composition[el_1]
            else:
                self.chem_composition[el_1] = 0.0

            m_2_n = 0
            if el_2 in self.chem_composition:
                m_2_n = self.chem_composition[el_2]
            else:
                self.chem_composition[el_2] = 0.0

            if el_2 == "C":
                m_i = mci
            elif el_2 == "N":
                m_i = mni

            mcarb = self.m_compound_calk(element="C", compounds=compounds, m_i=mci)
            mnitr = self.m_compound_calk(element="N", compounds=compounds, m_i=mni)

            marker_2 = el_1 + el_2
            k = self.eq_constant(composition=marker_2)

            f_1 = self.activity_coef(element=el_1, t=self.t, max_order=2)
            f_2 = self.activity_coef(element=el_2, t=self.t, max_order=2)
            f_compound = compounds_activity[el_1 + "_" + el_2]

            if mcarb == 0 and mnitr ==0:
                mcarb = epsilon
            else:
                mcarb = mcarb + epsilon

            elem_left =((m_elem_icompound) * (mr_1 + mr_2)) / ((mcarb + mnitr) * mr_1)
            elem_right = 10000.0 * ((m_2_n - m_i) / (100 - mcarb - mnitr)) * k * ((m_1_n - m_el_icarb - m_el_initr) /
                                                                              (100 - mcarb - mnitr)) * (
                      (f_1 * f_2) / (f_compound))

            if m_2_n == 0 or m_1_n == 0:
                return 0.0
            else:
                return (elem_left - elem_right) / (elem_left + epsilon)
        else:
            return 0


    def do_balance_equation(self, el_1="Ti", el_2="C", compounds={}, compounds_activity={}, initial_dif=10.0, minimal_el_cons=0.0, n_1=1, n_2=1):
        local_compounds = copy.deepcopy(compounds)

        rough_steps_n = 10.0
        i_1 = 0
        max_el_cons = self.max_element_consumption(el_1=el_1, el_2=el_2, n_1=float(n_1), n_2=float(n_2))
        d_i = (max_el_cons - minimal_el_cons) / float(rough_steps_n)
        current_dif = copy.deepcopy(initial_dif)
        dif_ctr = current_dif
        current_m_el_cons = minimal_el_cons
        current_m_el_cons_1 = minimal_el_cons
        current_m_el_cons_prev = minimal_el_cons

        current_dif = self.find_diff(el_1=el_1, el_2=el_2, compounds=local_compounds, compounds_activity=compounds_activity)
        dif_prev = 0.0 # copy.deepcopy(current_dif)

        iteration = 1
        while dif_ctr > self.convergention_ctiterial*10.0:
            i_1 = 0
            local_compounds[el_1 + "_" + el_2] = current_m_el_cons
            while i_1 < rough_steps_n + 1:

                current_dif = self.find_diff(el_1=el_1, el_2=el_2, compounds=local_compounds,
                                             compounds_activity=compounds_activity)

                if abs(current_dif) < self.convergention_ctiterial*10:
                    if abs(current_dif) <= dif_ctr:
                        current_m_el_cons_1 = current_m_el_cons
                        dif_ctr = current_dif

                if (current_dif > 0 and dif_prev < 0) or (current_dif < 0 and dif_prev > 0):
                    if (current_dif > 0 and dif_prev < 0 and d_i > 0) or (current_dif < 0 and dif_prev > 0 and d_i < 0):

                        current_m_el_cons_d = current_m_el_cons_prev
                        current_m_el_cons_t = current_m_el_cons

                        while abs(current_dif) > self.convergention_ctiterial:
                            current_m_el_cons = (current_m_el_cons_d + current_m_el_cons_t) / 2.0

                            local_compounds[el_1+"_"+el_2] = current_m_el_cons
                            current_dif = self.find_diff(el_1=el_1, el_2=el_2, compounds=local_compounds,
                                                         compounds_activity=compounds_activity)

                            if d_i > 0:
                                if current_dif < 0:
                                    current_m_el_cons_d = current_m_el_cons
                                elif current_dif > 0:
                                    current_m_el_cons_t = current_m_el_cons
                            else:
                                if current_dif > 0:
                                    current_m_el_cons_d = current_m_el_cons
                                elif current_dif < 0:
                                    current_m_el_cons_t = current_m_el_cons

                        if abs(current_dif) <= self.convergention_ctiterial*10:
                            if abs(current_dif) <= abs(dif_ctr):
                                current_m_el_cons_1 = current_m_el_cons
                                dif_ctr = current_dif

                    elif (current_dif > 0 and dif_prev < 0 and d_i < 0) or (current_dif < 0 and dif_prev > 0 and d_i > 0):
                        current_m_el_cons_t = current_m_el_cons_prev
                        current_m_el_cons_d = current_m_el_cons

                        while abs(current_dif) > self.convergention_ctiterial:
                            current_m_el_cons = (current_m_el_cons_d + current_m_el_cons_t) / 2.0
                            local_compounds[el_1 + "_" + el_2] = current_m_el_cons

                            current_dif = self.find_diff(el_1=el_1, el_2=el_2, compounds=local_compounds,
                                                         compounds_activity=compounds_activity)

                            if d_i > 0:
                                if current_dif < 0:
                                    current_m_el_cons_d = current_m_el_cons
                                elif current_dif > 0:
                                    current_m_el_cons_t = current_m_el_cons
                            else:
                                if current_dif > 0:
                                    current_m_el_cons_d = current_m_el_cons
                                elif current_dif < 0:
                                    current_m_el_cons_t = current_m_el_cons

                            #print(current_m_el_cons_d, current_m_el_cons_t)
                            #print(current_dif)

                        if abs(current_dif) <= self.convergention_ctiterial*10:
                            if abs(current_dif) <= abs(dif_ctr):
                                current_m_el_cons_1 = current_m_el_cons
                                dif_ctr = current_dif

                    i_1 = rough_steps_n

                dif_prev = current_dif
                current_m_el_cons_prev = current_m_el_cons
                if dif_ctr > self.convergention_ctiterial*10:
                    current_m_el_cons += d_i
                    if current_m_el_cons  < minimal_el_cons:
                        current_m_el_cons = minimal_el_cons
                local_compounds[el_1 + "_" + el_2] = current_m_el_cons

                i_1 += 1

            i_1 = rough_steps_n + 1

            if abs(dif_ctr) > self.convergention_ctiterial*10:

                iteration += 1
                rough_steps_n *= 10

                d_i = (max_el_cons - minimal_el_cons) / float(rough_steps_n)

                current_m_el_cons = minimal_el_cons
                current_m_el_cons_prev = current_m_el_cons
                local_compounds[el_1 + "_" + el_2] = current_m_el_cons
                i_1 = 0
            else:
                current_m_el_cons = current_m_el_cons_1

            if iteration > max_iterations:

                current_m_el_cons = max_el_cons
                current_dif = self.convergention_ctiterial/10.0
                dif_ctr = current_dif

        return (current_m_el_cons, current_dif, local_compounds)

    def moll_part_of_compound(self, el_1="Ti", el_2="C", mass=0.0, compounds={}):
        mci = self.m_i_calc(element="C", compounds=compounds)
        mni = self.m_i_calc(element="N", compounds=compounds)
        mcarb = self.m_compound_calk(element="C", compounds=compounds, m_i=mci)
        mnitr = self.m_compound_calk(element="N", compounds=compounds, m_i=mni)
        if mcarb <=0:
            mcarb = epsilon
        m_el_in_compound = mass
        mr_1 = 0
        if el_1 in self.chem.atom_mass and el_2 in self.chem.atom_mass:
            mr_1 = self.chem.atom_mass[el_1]
            mr_2 = self.chem.atom_mass[el_2]
            return (m_el_in_compound * (mr_1 + mr_2)) / ((mcarb + mnitr) * mr_1)
        return 0.0

    def percent_of_compound(self, moll_parts={}, compound="Ti_C"):
        if compound in moll_parts:
            x_el = moll_parts[compound]
            if compound.split("_")[0] in self.chem_composition and compound.split("_")[1] in self.chem_composition:
                mr_1 = self.chem_composition[compound.split("_")[0]]
                mr_2 = self.chem_composition[compound.split("_")[1]]
                factor = x_el * (mr_1 + mr_2)
                divider = 0.0
                for i in moll_parts:
                    x_current = moll_parts[i]
                    mr_1_current = 0
                    mr_2_current = 0
                    if i.split("_")[0] in self.chem_composition and i.split("_")[1] in self.chem_composition:
                        mr_1_current = self.chem_composition[i.split("_")[0]]
                        mr_2_current =self.chem_composition[i.split("_")[1]]
                    divider += x_current * (mr_1_current + mr_2_current)
                if divider > 0:
                     return 100.0 * (factor / divider)
        return float(0)

    def n_elements_calc(self, persents={}):
        compounds_n = {}
        for i in persents:
            if i.split("_")[0] in self.chem_composition and i.split("_")[1] in self.chem_composition:
                mr_1_current = self.chem_composition[i.split("_")[0]]
                mr_2_current = self.chem_composition[i.split("_")[1]]
                if (mr_1_current + mr_2_current) > 0:
                    compounds_n[i] = persents[i] / (mr_1_current + mr_2_current)
        elements_n = {}
        for i in self.chem_composition:
            for j in compounds_n:
                if i in j.split("_"):
                    if i not in elements_n:
                        elements_n[i] = compounds_n[j]
                    else:
                        elements_n[i] += compounds_n[j]
        return elements_n

    def psi_calc(self, element="Ti", x_elements={}, t=1000.0):
        if element in x_elements:
            sum = x_elements[element]
            for i in x_elements:
                if i != element:
                    if element in self.excange_energies:
                        if i in self.excange_energies[element]:
                            sum += x_elements[i] * math.exp(-self.excange_energies[element][i] / (ponamarenko_const * t))
            if sum <= 0:
                sum = epsilon
            return 1.0 / sum
        return 0.0

    def find_solvus_single(self, compound="Ti_C"):
        t_ptev = 0.0
        t_current = 1.0
        elements = compound.split("_")

        d_t = 10.0

        diff_prev = None
        diff = 100.0

        while diff**2 > self.convergention_ctiterial**2:

            g_full = self.calc_full_energy(compound=compound, t=t_current)
            lg_concentrations_sum = math.fsum([math.log10(self.chem_composition[i]) for i in elements])
            lg_activities_sum = math.fsum([math.log10(self.activity_coef(element=i, t=t_current)) for i in elements])

            left = g_full / (2.3 * R * t_current)
            right = lg_concentrations_sum + lg_activities_sum
            diff = (2 * (left - right) / (left + right))

            if (diff_prev is None) or (diff > 0 and diff_prev > 0) or (diff < 0 and diff_prev < 0):
                while (diff_prev is None) or (diff > 0 and diff_prev > 0) or (diff < 0 and diff_prev < 0):
                    t_ptev = t_current
                    t_current += d_t
                    diff_prev = diff
                    g_full = self.calc_full_energy(compound=compound, t=t_current)
                    lg_concentrations_sum = math.fsum([math.log10(self.chem_composition[i]) for i in elements])
                    lg_activities_sum = math.fsum([math.log10(self.activity_coef(element=i, t=t_current)) for i in elements])

                    left = g_full/(2.3 * R * t_current)
                    right = lg_concentrations_sum + lg_activities_sum
                    diff = (2*(left - right)/(left + right))
                    print(t_current, diff)

            t_current = t_ptev - diff_prev * ((t_ptev - t_current)/(diff_prev - diff))
            print(t_current, diff)

        return  t_current

    def find_solvus_carbonitride(self, element="Ti"):
        t_current = 1873.0

        diff_t = None
        while diff_t is None or diff_t**2 > self.convergention_ctiterial**2:
            lg_activity_element = math.log10(self.activity_coef(element=element, t=t_current))
            lg_activity_c = math.log10(self.activity_coef(element="C", t=t_current))
            lg_activity_n = math.log10(self.activity_coef(element="N", t=t_current))

            s_nitr = 0.0
            s_carb = 0.0
            q_nitr = 0.0
            q_carb = 0.0
            if "eq_constant_params" in self.td_params:
                if element + "N" in self.td_params["eq_constant_params"]:
                    q_nitr = self.td_params["eq_constant_params"][element + "N"]["Q"]
                    s_nitr = self.td_params["eq_constant_params"][element + "N"]["S"]
                if element + "C" in self.td_params["eq_constant_params"]:
                    q_carb = self.td_params["eq_constant_params"][element + "C"]["Q"]
                    s_carb = self.td_params["eq_constant_params"][element + "C"]["S"]

            nitr_fctor = (math.log10(self.chem_composition[element]*self.chem_composition["N"]) +
                         lg_activity_element + lg_activity_n - s_nitr)

            carb_fctor = (math.log10(self.chem_composition[element] * self.chem_composition["C"])
                         + lg_activity_element + lg_activity_c - s_carb)

            x_carb = 0.0
            diff = None
            min_eps = 0.1
            counter = 1.0
            while diff is None or diff**2 > self.convergention_ctiterial**2:
                x_carb_new = math.exp(2.3 * (carb_fctor - (q_carb/q_nitr) * nitr_fctor)) * pow((abs(1.0 - x_carb)), (q_carb/q_nitr))
                diff = 2* (x_carb - x_carb_new) / (x_carb_new + x_carb + epsilon)
                x_carb = x_carb_new

                if x_carb > 1.0:
                    x_carb = 0.9999
                    break

                # print(diff, x_carb)

            x_nitr = 1.0 - x_carb

            t_current_new = q_nitr/(math.log10(x_nitr) - nitr_fctor)
            diff_t = 2* (t_current - t_current_new) / (t_current_new + t_current + epsilon)
            t_current = t_current_new

        return t_current

    def solve(self):
        initial_carbides_and_nitrides_concentrations = {
            "Ti_C": 0,
            "Nb_C": 0,
            "V_C": 0,
            "Al_N": 0,
            "Zr_C": 0,
            "Ti_N": 0,
            "Nb_N": 0,
            "V_N": 0,
            "Zr_N": 0
        }

        initial_activity_for_compounds = {
            "Ti_C": 1,
            "Nb_C": 1,
            "V_C": 1,
            "Al_N": 1,
            "Zr_C": 1,
            "Ti_N": 1,
            "Nb_N": 1,
            "V_N": 1,
            "Zr_N": 1
        }

        dif_compounds ={
            "Ti_C": 100,
            "Nb_C": 100,
            "V_C": 100,
            "Al_N": 100,
            "Zr_C": 100,
            "Ti_N": 100,
            "Nb_N": 100,
            "V_N": 100,
            "Zr_N": 100
        }

        percents_of_compounds_final = {}

        current_carbides_and_nitrides_concentrations = copy.deepcopy(initial_carbides_and_nitrides_concentrations)

        global_iteration_criteria = [100.0 for i in initial_carbides_and_nitrides_concentrations]
        global_iterator = 0
        while max(global_iteration_criteria) > self.convergention_ctiterial*100 and global_iterator < max_iterations:
            global_iterator += 1
            for i in initial_carbides_and_nitrides_concentrations:
                initial_carbides_and_nitrides_concentrations[i] = 0

            for i in dif_compounds:
                dif_compounds[i] = 100

            old_carbides_and_nitrides_concentrations = copy.deepcopy(current_carbides_and_nitrides_concentrations)
            current_carbides_and_nitrides_concentrations = copy.deepcopy(initial_carbides_and_nitrides_concentrations)
            current_activity_for_compounds = copy.deepcopy(initial_activity_for_compounds)
            dif_compounds_ctr = copy.deepcopy(dif_compounds)
            square_diff = [abs(dif_compounds_ctr[i]) for i in dif_compounds_ctr]
            max_key_dif = ""
            i_tmp = 1
            while max(square_diff) > self.convergention_ctiterial*1000:

                for i in initial_activity_for_compounds:
                    el_1, el_2 = i.split("_")
                    if i in dif_compounds and i != max_key_dif:
                        dif_compounds[i] = self.find_diff(el_1=el_1, el_2=el_2, compounds=current_carbides_and_nitrides_concentrations,
                                                          compounds_activity=current_activity_for_compounds)
                i_tmp += 1
                dif_compounds_ctr = copy.deepcopy(dif_compounds)
                dif_compounds_square = copy.deepcopy(dif_compounds)
                for i in dif_compounds_square:
                    dif_compounds_square[i] = dif_compounds[i] ** 2

                max_key_dif, max_dif = max(dif_compounds_square.items(), key=operator.itemgetter(1))

                el_1_current, el_2_current = max_key_dif.split("_")
                tmp_new_conc, tmp_diff, tmp_compounds = self.do_balance_equation(el_1=el_1_current, el_2=el_2_current, compounds=current_carbides_and_nitrides_concentrations,
                                    compounds_activity=initial_activity_for_compounds)

                current_carbides_and_nitrides_concentrations[max_key_dif] = tmp_new_conc
                dif_compounds[max_key_dif] = tmp_diff
                square_diff = [abs(dif_compounds[i]) for i in dif_compounds]

            x_compounsds = {}

            for i in current_carbides_and_nitrides_concentrations:
                el_1 = i.split("_")[0]
                el_2 = i.split("_")[1]
                x_compounsds[i] = self.moll_part_of_compound(el_1=el_1, el_2=el_2,
                                                            mass=current_carbides_and_nitrides_concentrations[i],
                                                            compounds=current_carbides_and_nitrides_concentrations)

            percents_of_compounds = {}
            for i in x_compounsds:
                percents_of_compounds[i] = self.percent_of_compound(moll_parts=x_compounsds, compound=i)

            elements_n = self.n_elements_calc(persents=percents_of_compounds)
            summ = 0
            for i in elements_n:
                summ += elements_n[i]

            x_elements = {}
            for i in elements_n:
                if summ != 0:
                    x_elements[i] = elements_n[i]/summ
                else:
                    x_elements[i] = 0.0

            psi_values = {}
            for i in x_elements:
                psi_values[i] = self.psi_calc(element=i, x_elements=x_elements, t=self.t)

            for i in initial_activity_for_compounds:
                el_1 = i.split("_")[0]
                el_2 = i.split("_")[1]
                if el_1 in psi_values and el_2 in psi_values:
                    initial_activity_for_compounds[i] = psi_values[el_1] * psi_values[el_2]
                    # print(i, initial_activity_for_compounds[i], el_1, el_2, psi_values[el_1], psi_values[el_2])

            global_iteration_criteria_new = []

            for i in old_carbides_and_nitrides_concentrations:
                if i in old_carbides_and_nitrides_concentrations and i in current_carbides_and_nitrides_concentrations:
                    global_iteration_criteria_new.append(((old_carbides_and_nitrides_concentrations[i] - current_carbides_and_nitrides_concentrations[i])/(
                    old_carbides_and_nitrides_concentrations[i] + epsilon))**2)
            global_iteration_criteria = global_iteration_criteria_new

            percents_of_compounds_final = copy.deepcopy(percents_of_compounds)

        mci = self.m_i_calc(element="C", compounds=current_carbides_and_nitrides_concentrations)
        mni = self.m_i_calc(element="N", compounds=current_carbides_and_nitrides_concentrations)
        mcarb = self.m_compound_calk(element="C", compounds=current_carbides_and_nitrides_concentrations, m_i=mci)
        mnitr = self.m_compound_calk(element="N", compounds=current_carbides_and_nitrides_concentrations, m_i=mni)

        elements_in_solution = self.calc_solution_composition(current_compuunds_concentrations=current_carbides_and_nitrides_concentrations)

        ounput = {
            "elements_consumption": current_carbides_and_nitrides_concentrations,
            "percents_of_compounds": percents_of_compounds_final,
            "carbides_mass": mcarb,
            "nitrides_mass": mnitr,
            "elements_in_solution": elements_in_solution
        }

        return ounput

    def calc_solution_composition(self, current_compuunds_concentrations={}):
        mci = self.m_i_calc(element="C", compounds=current_compuunds_concentrations)
        mni = self.m_i_calc(element="N", compounds=current_compuunds_concentrations)
        mcarb = self.m_compound_calk(element="C", compounds=current_compuunds_concentrations, m_i=mci)
        mnitr = self.m_compound_calk(element="N", compounds=current_compuunds_concentrations, m_i=mni)

        c_percent = 100 * (self.chem_composition["C"] - mci) / (100 - mcarb - mnitr)
        n_percent = 100 * (self.chem_composition["N"] - mni) / (100 - mcarb - mnitr)

        elements_in_solution = {}

        for i in self.chem_composition:
            pr_tmp = self.chem_composition[i]
            consumption = 0
            for j in current_compuunds_concentrations:
                if i in j.split("_"):
                    consumption += current_compuunds_concentrations[j]
            if consumption > 0:
                elements_in_solution[i] = 100 * (pr_tmp - consumption) / (100 - mcarb - mnitr)
            else:
                elements_in_solution[i] = pr_tmp

        elements_in_solution["N"] = n_percent
        elements_in_solution["C"] = c_percent

        return elements_in_solution

'''
print(cbn.check_carbonitride_formers())
print(cbn.calculate_excange_energies())
for i in chemComposition:
    print(i, cbn.activity_coef(element=i, t=1300))
print()
print("Ti", cbn.max_element_consumption(el_1="Ti"))
print("V", cbn.max_element_consumption(el_1="V"))
print("Nb", cbn.max_element_consumption(el_1="Nb"))
'''
'''
initial_carbides_and_nitrides_concentrations = {
    "Ti_C": 0,
    "Nb_C": 0,
    "V_C": 0,
    "Al_N": 0,
    "Zr_C": 0,
    "Ti_N": 0,
    "Nb_N": 0,
    "V_N": 0,
    "Zr_N": 0
}

initial_activity_for_compounds = {
    "Ti_C": 1,
    "Nb_C": 1,
    "V_C": 1,
    "Al_N": 1,
    "Zr_C": 1,
    "Ti_N": 1,
    "Nb_N": 1,
    "V_N": 1,
    "Zr_N": 1
}

dif_compounds ={
    "Ti_C": 100,
    "Nb_C": 100,
    "V_C": 100,
    "Al_N": 100,
    "Zr_C": 100,
    "Ti_N": 100,
    "Nb_N": 100,
    "V_N": 100,
    "Zr_N": 100
}

current_carbides_and_nitrides_concentrations = copy.deepcopy(initial_carbides_and_nitrides_concentrations)

#print(cbn.do_balance_equation(el_1="Ti", el_2="C", compounds=current_carbides_and_nitrides_concentrations,
#                        compounds_activity=initial_activity_for_compounds))

for t in range(900, 1600, 40):
    cbn = CarbonitrideThermodynamicsInSteel(chem_composition=chemComposition, td_params=td_params, t=t)
    print(t-273, cbn.solve())
'''

#cbn = CarbonitrideThermodynamicsInSteel(chem_composition=chemComposition, td_params=td_params, t=1200)
'''
print(cbn.calc_solute_energy(element="C", t=1200.0))
print(cbn.calc_solute_energy(element="Ti", t=1000.0))
print(cbn.calc_solute_energy(element="Nb", t=1300.0))

print(cbn.calc_full_energy(compound="Nb_C", t=390.0))
print(cbn.calc_full_energy(compound="Nb_C", t=800.0))
print(cbn.calc_full_energy(compound="Nb_C", t=1100.0))


print(cbn.find_solvus_carbonitride(element="Ti"))
print(cbn.find_solvus_carbonitride(element="Nb"))
print(cbn.find_solvus_carbonitride(element="V"))
print(cbn.find_solvus_carbonitride(element="Al"))
print(cbn.find_solvus_carbonitride(element="Zr"))

'''

