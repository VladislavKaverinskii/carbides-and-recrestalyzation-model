# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
import flow_stress
import rec_press
import thermodytamic
import copy
import functools

R = 8.3144598
N_A = 6.0221409e+23
K_B = 1.3806485279e-23

disl_dens_test = 1892308866164430.5


class Nucleation:
    def __init__(self, parameters, composition, init_t=1000.0, base="Fe", lattice_type="FCC"):
        self.parameters = parameters
        self.composition = composition
        self.current_t = init_t

        self.disl_handler = rec_press.Dislocations(parameters=self.parameters, composition=composition)

        self.chem = rec_press.Chemical(composition, base=base, lattice_type=lattice_type)

        self.thermo_solver = None
        if base == "Fe":
            self.thermo_solver = thermodytamic.CarbonitrideThermodynamicsInSteel(chem_composition=composition,
                                                                                 td_params = self.parameters["td_params"],
                                                                                 t = init_t)

        self.equilibrium_data = self.calculate_carbonitrides_equilibrium()

        self.compounds_compositions = {}
        if self.equilibrium_data is not None:
            self.compounds_compositions = self.calc_carbonitrides_stoichiometry()

    def reset_temperatute(self, new_t):
        if self.thermo_solver is not None:
            self.thermo_solver = thermodytamic.CarbonitrideThermodynamicsInSteel(chem_composition=self.composition,
                                                                                 td_params=self.parameters["td_params"],
                                                                                 t=new_t)

    def calculate_carbonitrides_equilibrium(self):
        if self.thermo_solver is not None:
            return self.thermo_solver.solve()
        return None

    @functools.lru_cache(maxsize=128, typed=False)
    def calc_carbonitrides_stoichiometry(self):
        if self.equilibrium_data is not None:
            if "percents_of_compounds" in self.equilibrium_data:
                carbonitrides_masses = {}
                for i in self.equilibrium_data["percents_of_compounds"]:
                    first_element = i.split("_")[0]
                    if len(i.split("_")) > 1:
                        if i.split("_")[1] == "C" or i.split("_")[1] == "N":
                            if first_element+"_C"+"_N" not in carbonitrides_masses:
                                carbonitrides_masses[first_element+"_C"+"_N"] = self.equilibrium_data["percents_of_compounds"][i]
                            else:
                                carbonitrides_masses[first_element + "_C" + "_N"] += self.equilibrium_data["percents_of_compounds"][i]
                carbonitrides_compositions = {}
                for i in self.equilibrium_data["percents_of_compounds"]:
                    elements = i.split("_")
                    for j in carbonitrides_masses:
                        carbonitride_elements = j.split("_")
                        if j not in carbonitrides_compositions:
                            carbonitrides_compositions[j] = {}
                            carbonitrides_compositions[j][carbonitride_elements[0]] = 1.0
                            if elements[0] == carbonitride_elements[0]:
                                if "N" in elements:
                                    carbonitrides_compositions[j]["N"] = (self.equilibrium_data["percents_of_compounds"][i] /
                                                                      carbonitrides_masses[j])
                                if "C" in elements:
                                    carbonitrides_compositions[j]["C"] = (self.equilibrium_data["percents_of_compounds"][i] /
                                                                      carbonitrides_masses[j])
                        else:
                            if elements[0] == carbonitride_elements[0]:
                                if "N" in elements:
                                    carbonitrides_compositions[j]["N"] = (self.equilibrium_data["percents_of_compounds"][i] /
                                                                      carbonitrides_masses[j])
                                if "C" in elements:
                                    carbonitrides_compositions[j]["C"] = (self.equilibrium_data["percents_of_compounds"][i] /
                                                                      carbonitrides_masses[j])
                for i in carbonitrides_compositions:
                    if "C" not in carbonitrides_compositions[i]:
                        carbonitrides_compositions[i]["C"] = 0
                    if  "N" not in carbonitrides_compositions[i]:
                        carbonitrides_compositions[i]["N"] = 0
                return carbonitrides_compositions
        return {}

    @functools.lru_cache(maxsize=128, typed=False)
    def conpounds_molar_volumes(self):
        molar_volumes = {}
        compounds_compositions = {}
        if self.compounds_compositions is not None:
             compounds_compositions = self.compounds_compositions
        elif "compounds_compositions" in self.parameters:
            compounds_compositions = self.parameters["compounds_compositions"]
        else:
            return compounds_compositions
        compounds_dens = self.parameters["compounds_dens"]
        # it assumes that density is a constant and not depends on composition

        for i in compounds_compositions:
            m_r = 0
            list_of_semicompounds = []
            dens = 0
            parts = {}
            sum_part = 0
            for j in compounds_compositions[i]:
                m_r += compounds_compositions[i][j] * self.chem.atom_mass[j] / 1000.0
                for k in compounds_dens:
                    if j in k.split("_"):
                        for j_2 in compounds_compositions[i]:
                            if j_2 != j and j_2 in k.split("_"):
                                if k not in list_of_semicompounds:
                                    list_of_semicompounds.append(k)
                                    part = min(compounds_compositions[i][j], compounds_compositions[i][j_2])
                                    parts[k] = part
                                    sum_part += part

            for p in parts:
                dens += (parts[p] / sum_part) * compounds_dens[p]

            molar_volumes[i] = m_r / dens
        return molar_volumes

    @functools.lru_cache(maxsize=128, typed=False)
    def calc_alpha_n(self, element="Nb"):
        el_c = 0
        el_n = 0
        c_el = thermodytamic.epsilon
        n_el = thermodytamic.epsilon
        if "td_params" in self.parameters:
            if "vagner_interactions" in self.parameters["td_params"]:
                interaction_parameters = self.parameters["td_params"]["vagner_interactions"]
                if element in interaction_parameters:
                    if "e" in interaction_parameters[element]:
                        if "C" in interaction_parameters[element]["e"]:
                            el_c = interaction_parameters[element]["e"]["C"]["s"]
                        if "N" in interaction_parameters[element]["e"]:
                            el_n = interaction_parameters[element]["e"]["N"]["s"]
                if "C" in interaction_parameters:
                    if "e" in interaction_parameters["C"]:
                        if element in interaction_parameters["C"]["e"]:
                            c_el = interaction_parameters["C"]["e"][element]["s"]
                if "N" in interaction_parameters:
                    if "e" in interaction_parameters["N"]:
                        if element in interaction_parameters["N"]["e"]:
                            n_el = interaction_parameters["N"]["e"][element]["s"]

        return abs((el_c + el_n)/(c_el + n_el))

    def form_list_of_alpha_n(self):
        list_of_alpha_n = {}
        if self.compounds_compositions is not None:
            for i in self.compounds_compositions:
                list_of_alpha_n[i] = self.calc_alpha_n(i.split("_")[0])
        return list_of_alpha_n

    def calc_td_forces(self, current_solution={}, t=1000.0):
        v_m_list = self.conpounds_molar_volumes()
        if self.equilibrium_data is not None:
            if "elements_in_solution" in self.equilibrium_data:
                eq_elements_in_solution = self.equilibrium_data["elements_in_solution"]
                eq_x_element_list, eq_mol_of_solution = self.calc_molar_concentrations(eq_elements_in_solution)
                alpha_n_list = self.form_list_of_alpha_n()
                eq_x_c_effective_list = {}
                for i in self.compounds_compositions:
                    c_i = self.equilibrium_data["elements_in_solution"]["C"]
                    n_i = self.equilibrium_data["elements_in_solution"]["N"]
                    m_c_eff = c_i + alpha_n_list[i] * n_i
                    moll_c_eff = m_c_eff / self.chem.atom_mass["C"]
                    x_c_eff = moll_c_eff / eq_mol_of_solution
                    eq_x_c_effective_list[i] = x_c_eff
                current_elements_in_solution = current_solution
                current_x_element_list, current_mol_of_solution = self.calc_molar_concentrations(current_elements_in_solution)
                current_x_c_effective_list = {}
                for i in self.compounds_compositions:
                    c_i = current_solution["C"]
                    n_i = current_solution["N"]
                    m_c_eff = c_i + alpha_n_list[i] * n_i
                    moll_c_eff = m_c_eff / self.chem.atom_mass["C"]
                    x_c_eff = moll_c_eff / current_mol_of_solution
                    current_x_c_effective_list[i] = x_c_eff
                g_p_list = {}
                for i in v_m_list:
                    element = i.split("_")[0]
                    current_x_element = current_x_element_list[element]
                    current_x_c = current_x_c_effective_list[i]
                    eq_x_element = eq_x_element_list[element]
                    eq_x_c = eq_x_c_effective_list[i]
                    '''
                    print("--v_m_list--")
                    print(i)
                    print("v_m_list[i]", v_m_list[i])
                    print("eq_x_element", eq_x_element)
                    print("eq_x_c", eq_x_c)
                    print(current_x_element * current_x_c)
                    print(eq_x_element * eq_x_c)
                    print((current_x_element * current_x_c)/(eq_x_element * eq_x_c))
                    print(math.log((current_x_element * current_x_c)/(eq_x_element * eq_x_c)))
                    print("-----------------")
                    '''
                    #print(t)
                    g_p_list[i] = - ((R * t) / (v_m_list[i] * 7.0)) * math.log((current_x_element * current_x_c)/
                                                                       (eq_x_element * eq_x_c))


                return g_p_list

    def calc_molar_concentrations(self, elements_concentrations):
        elements_in_solution = elements_concentrations
        mols_of_components = {}
        for i in elements_in_solution:
            mols_of_components[i] = elements_in_solution[i]/self.chem.atom_mass[i]
        mol_of_solution = math.fsum([elements_in_solution[i]/self.chem.atom_mass[i]
                              for i in elements_in_solution])
        molar_concentrations = {}
        for i in mols_of_components:
            molar_concentrations[i] = mols_of_components[i] / mol_of_solution
        return molar_concentrations, mol_of_solution

    def n_tau(self, x_t=0.0, D_0=1.0, e=0.2, v=1.0, T=1000.0, element="Nb", T_0=1000.0, tau=0.0, n_total=0):
        c = self.composition["C"]
        disl_dens = self.disl_handler.austenite_disl_dens_t(x_t=x_t, D_0=D_0, e=e, v=v, T=T, C=c, tau=tau, n_total=n_total)
        n_n = self.n_n(disl_dens=disl_dens, n=n_total)
        d = self.diffusion_coef(t=T_0, element=element)
        a = self.disl_handler.austenite_a(T=T, C=c)

    @functools.lru_cache(maxsize=128, typed=False)
    def n_n(self, disl_dens=1.0, n=0):
        return 0.5 * pow(disl_dens, 1.5) - n

    def beta(self, compound="Nb_C_N", t=1000.0, current_solution={}):
        a_gamma = self.disl_handler.austenite_a(T=t, C=current_solution["C"])
        #a_gamma = self.gamma(t=t)
        #print("t, beta", t)
        r_c = self.r_c(t=t, compound=compound, current_solution=current_solution)
        element=compound.split("_")[0]
        if element in current_solution:
            x_element=current_solution[element]
            d = self.diffusion_coef(element=element, t=t)
            return (4.0*math.pi*(r_c**2)*d*x_element)/(a_gamma**4)
        else:
            return 0.0

    @functools.lru_cache(maxsize=128, typed=False)
    def gamma_p_g(self, element="Nb"):
        compound = element + "_C_N"
        t_s = self.thermo_solver.find_solvus_carbonitride(element=element)
        if "surface_energy" in self.parameters:
            if compound in self.parameters["surface_energy"]:
                intercept = self.parameters["surface_energy"][compound]["intercept"]
                pre_exp = self.parameters["surface_energy"][compound]["pre_exp"]
                exp_factor = self.parameters["surface_energy"][compound]["exp_factor"]
                return intercept + pre_exp * math.exp(exp_factor * t_s)
        return 0.09272 + 5.84e25 * math.exp(-0.04416 * t_s)

    def r_c(self, t=1000.0, compound="Nb_C_N", current_solution={}):
        element = compound.split("_")[0]
        g = self.gamma_p_g(element=element)
        #print("t, r_c", t)
        g_p_list = self.calc_td_forces(current_solution=current_solution, t=t)
        g_p = thermodytamic.epsilon
        if compound in g_p_list:
            g_p = g_p_list[compound]

        #print("g", g)
        #print("g_p", g_p)
        return -(2.0 * g) / g_p

    def g_c(self, t=1000.0, compound="Nb_C_N", current_solution={}):
        #rint("t, g_c", t)
        r = self.r_c(t=t, compound=compound, current_solution=current_solution)
        element = compound.split("_")[0]
        g = self.gamma_p_g(element=element)
        #print("--g_c--")
        #print(element)
        #print("r", r)
        #print("g", g)
        #print("----------")
        return (4.0 / 3.0) * math.pi * r * r * g

    def zeldovich_factor(self, t=1000.0, compound="Nb_C_N", current_solution={}):
        g = self.g_c(t=t, compound=compound, current_solution=current_solution)
        n = self.n_c(t=t, compound=compound, current_solution=current_solution)
        return pow((g / (3.0 * math.pi * K_B * t * (n**2))), 0.5)

    @functools.lru_cache(maxsize=128, typed=False)
    def diffusion_coef(self, element="Nb", t=1000.0):
        if element in self.parameters["diffusion_parameters"]:
            a = self.parameters["diffusion_parameters"][element]["a"]
            q = self.parameters["diffusion_parameters"][element]["q"]
            return a*math.exp(q/(R*t))
        else:
            return 0.0

    def n_c(self, t=1000.0, compound="Nb_C_N", current_solution={}):
        v_m = self.conpounds_molar_volumes()[compound]
        v_a = v_m / N_A
        #print("t, n_c", t)
        r = self.r_c(t=t, compound=compound, current_solution=current_solution)
        return (4.0 / 3.0) * math.pi * (r ** 3) / v_a

    def nucleation_rates(self, current_solution={}, t=1000.0, disl_dens=1.0, n=0):
        #print("t, nucleation_rates", t)
        g_of_compounds = self.calc_td_forces(current_solution=current_solution, t=t)

        g_c_list = {}
        for i in g_of_compounds:
            element = i.split("_")[0]
            t_s = self.thermo_solver.find_solvus_carbonitride(element=element)
            if t < t_s:
                g_c_list[i] = self.g_c(t=t, compound=i, current_solution=current_solution)

        n_n = self.n_n(disl_dens=disl_dens, n=n)
        #print(n_n)
        output = {}
        for i in g_c_list:
            z_factor = self.zeldovich_factor(t=t, compound=i, current_solution=current_solution)
            beta = self.beta(t=t, compound=i, current_solution=current_solution)
            '''
            print('-------')
            print(i)
            print('n_n', n_n)
            print('z_factor', z_factor)
            print('beta', beta)
            print('g_c_list[i]', g_c_list[i])
            print('math.exp(-g_c_list[i]/(R*t))', math.exp(-g_c_list[i]/(K_B*t)))
            print('n_rate', n_n * z_factor * beta * math.exp(-g_c_list[i]/(K_B*t)))
            print('-------')
            '''
            output[i] = n_n * z_factor * beta * math.exp(-g_c_list[i]/(K_B*t))
        return output


class Growth:
    def __init__(self, parameters, composition, init_t=1350.0, base="Fe", lattice_type="FCC"):
        self.nucl_handler = Nucleation(parameters, composition, init_t=init_t)
        self.parameters = parameters
        self.composition = composition
        self.current_t = init_t
        self.chem = rec_press.Chemical(composition, base=base, lattice_type=lattice_type)
        self.equilibrium_data = self.nucl_handler.equilibrium_data
        self.v_m_list = self.nucl_handler.conpounds_molar_volumes()
        self.a_c = 0.03

        for i in self.equilibrium_data:
            print(i, self.equilibrium_data[i])

        print(self.nucl_handler.compounds_compositions)

    @functools.lru_cache(maxsize=128, typed=False)
    def element_dens(self, element="Fe", t=1000.0):
        if "elements_dens" in self.parameters:
            if element in self.parameters["elements_dens"]:
                a = 0
                b = 0
                c = 0
                if "a" in self.parameters["elements_dens"][element]:
                    a = self.parameters["elements_dens"][element]["a"]
                if "b" in self.parameters["elements_dens"][element]:
                    b = self.parameters["elements_dens"][element]["b"]
                if "c" in self.parameters["elements_dens"][element]:
                    c = self.parameters["elements_dens"][element]["c"]
                return (a * (t**2)) + (b * t) + c
        return 0.0

    def solution_dens(self, t=1000.0, composition={}):
        d = 0.0
        for element in composition:
            d += self.element_dens(element=element, t=t) * composition[element] / 100.0
        return d

    def austenite_molar_volume(self, t=1000.0, composition={}):
        return (self.chem.MR / 1000.0) / self.solution_dens(t=t, composition=composition)

    def alpha_p(self, compound="Nb_C_N", composition={}, t=1000.0):
        if compound in self.v_m_list:
            return self.austenite_molar_volume(t=t, composition=composition) / self.v_m_list[compound]
        return 0.0

    def r_0_back(self, composition="Nb_C_N", t=1000.0):
        # do not use!!!
        if composition in self.v_m_list:
            element = composition.split("_")[0]
            gamma = self.nucl_handler.gamma_p_g(element=element)
            v = self.v_m_list[composition]
            return (2.0 * gamma * v) / (R * t)
        return 0.0

    def r_c_eff(self, t=1000.0, compound="Nb_C_N", current_solution={}):
        #print("t, r_c_eff", t)
        return (1.0 + self.a_c) * self.nucl_handler.r_c(t=t, compound=compound, current_solution=current_solution)

    def austenite_mr(self, solution_composition):
        handler = rec_press.Chemical(solution_composition)
        return handler.MR

    def eq_x_concentration(self, element="Nb"):
        if "elements_in_solution" in self.equilibrium_data:
            if element in self.equilibrium_data["elements_in_solution"]:
                mr = self.austenite_mr(self.equilibrium_data["elements_in_solution"])
                return (self.chem.atom_mass[element] * self.equilibrium_data["elements_in_solution"][element]) / mr
        return 0.0

    def actual_x_concentration(self, element="Nb", current_solution={}):
        if element in current_solution:
            mr = self.austenite_mr(current_solution)
            return (self.chem.atom_mass[element] * current_solution[element]) / mr
        return 0.0

    def x_p_element(self, element="Nb"):
        # part of an element in carbonitride of the element
        composition_name = element + "_C_N"
        el_part = self.nucl_handler.compounds_compositions[composition_name][element]
        c_part = self.nucl_handler.compounds_compositions[composition_name]["C"]
        n_part = self.nucl_handler.compounds_compositions[composition_name]["N"]
        return el_part / (el_part + c_part + n_part)

    @functools.lru_cache(maxsize=128, typed=False)
    def gamma_p_gamma(self, compound="Nb_C_N", t_sovus=1000.0):
        if "surface_energy" in self.parameters:
            if compound in self.parameters["surface_energy"]:
                intercept = self.parameters["surface_energy"][compound]["intercept"]
                pre_exp = self.parameters["surface_energy"][compound]["pre_exp"]
                exp_factor = self.parameters["surface_energy"][compound]["exp_factor"]
                return intercept + pre_exp * math.exp(exp_factor * t_sovus)
        return 0.09272 + 5.84e25 * math.exp(-0.04416 * t_sovus)

    @functools.lru_cache(maxsize=128, typed=False)
    def r_0(self, element="Nb", t=1350.0):
        composition_name = element + "_C_N"
        if composition_name in self.v_m_list:
            t_s = self.nucl_handler.thermo_solver.find_solvus_carbonitride(element=element)
            #if element=="V":
            #    print(t, t_s)
            g = self.gamma_p_gamma(compound=composition_name, t_sovus=t_s)
            v_mol = self.v_m_list[composition_name]
            return (2.0 * g * v_mol) / (R * t)
        return 0.0

    @functools.lru_cache(maxsize=128, typed=False)
    def disl_diff_coef(self, elemnt="Nb", t=1350.0):
        if "dislocations_diffusion" in self.parameters:
            if elemnt in self.parameters["dislocations_diffusion"]:
                a = self.parameters["dislocations_diffusion"][elemnt]["a"]
                q = self.parameters["dislocations_diffusion"][elemnt]["q"]
                return a * math.exp(q/(R*t))
        return 0.0

    @functools.lru_cache(maxsize=128, typed=False)
    def diff_coef(self, elemnt="Nb", t=1350.0):
        if "diffusion_parameters" in self.parameters:
            if elemnt in self.parameters["diffusion_parameters"]:
                a = self.parameters["diffusion_parameters"][elemnt]["a"]
                q = self.parameters["diffusion_parameters"][elemnt]["q"]
                return a * math.exp(q / (R * t))
        return 0.0

    def diff_coef_eff(self, r_current=0.0, elemnt="Nb", t=1350.0, disl_dens=0.0, current_solution={}):
        c = current_solution["C"]
        a = (self.nucl_handler.disl_handler.austenite_a(T=t, C=c),)
        b = self.nucl_handler.disl_handler.burgers(params=a)
        r_d = 3.0*b[0]
        d = self.diff_coef(elemnt=elemnt, t=t)
        d_pipe = self.disl_diff_coef(elemnt=elemnt, t=t)

        if r_current > 0:
            r_coef = math.exp(-0.00000000883569754/r_current)
        else:
            r_coef = 0.0

        if r_current < 5e9:
            r_coef = 0.0

        #print("r_d", r_d)
        #print("disl_dens", disl_dens)
        #print("d_pipe", d_pipe)
        #print("d", d)
        #print("math.pi * (r_d**2.0) * disl_dens", math.pi * (r_d**2.0) * disl_dens)
        return r_coef * math.pi * (r_d**2.0) * disl_dens * d_pipe + (1 - r_coef * math.pi * (r_d**2.0) * disl_dens) * d

    def growth_rate(self, element="Nb", d_tau=0.001, t=1350.0, disl_dens=0.0, n_current=0.0,
                    r_current=1e-10, v_nucl=0.0, current_solution={}):
        d = self.diff_coef_eff(elemnt=element, r_current=r_current, t=t, disl_dens=disl_dens,
                               current_solution=current_solution)
        r_0 = self.r_0(element=element, t=t)
        x_elem = self.actual_x_concentration(element=element, current_solution=current_solution)
        x_eq = self.eq_x_concentration(element=element)
        x_p = self.x_p_element(element=element)
        compound = element + "_C_N"
        #print("t, growth_rate", t)
        r_c = self.r_c_eff(t=t, compound=compound, current_solution=current_solution)
        '''
        print(compound)
        print("r_c ", r_c)
        print("r_0 ", r_0)
        print("r_current", r_current)
        print("d", d)
        print("x_elem", x_elem)
        print("x_eq", x_eq)
        print("x_p", x_p)
        print("v_nucl", v_nucl)
        print("n_current", n_current)
        '''
        a = self.alpha_p(compound=compound, composition=current_solution, t=t)

        #print("a", a)

        try:
            part_1 = d / r_current
        except:
            part_1 = 0

        part_2 = (x_elem - (x_eq * math.exp(r_0 / r_current))) / ((a * x_p) - (x_eq * math.exp(r_0 / r_current)))


        part_3 = (v_nucl / n_current) * (r_0 - r_current)
        # part_3 = (r_0 * v_nucl - r_current * n_current) / n_current

        if part_1 < 0:
            part_1 = 0

        #if part_2 < 0:
        #    part_2 = (x_elem - (x_eq * math.exp(-r_0 / r_current))) / ((a * x_p) - (x_eq * math.exp(-r_0 / r_current)))

        #print(compound)
        '''
        print('part_1', part_1)
        print('part_2', part_2)
        print('part_3', part_3)
        print("g_rate", (part_1 * part_2) + part_3)
        '''


        return (part_1 * part_2) + part_3

    def new_radius(self, element="Nb", d_tau=0.001, t=1350.0, disl_dens=0.0, d_n=0.0, n_current=0.0,
                    r_current=1e-10, v_nucl=0.0, current_solution={}):
        d = self.diff_coef_eff(elemnt=element, r_current=r_current, t=t, disl_dens=disl_dens,
                               current_solution=current_solution)
        r_0 = self.r_0(element=element, t=t)
        x_elem = self.actual_x_concentration(element=element, current_solution=current_solution)
        x_eq = self.eq_x_concentration(element=element)
        x_p = self.x_p_element(element=element)
        compound = element + "_C_N"
        #print("t, new_radius", t)
        r_c = self.r_c_eff(t=t, compound=compound, current_solution=current_solution)
        '''
        print(compound)
        print("r_c ", r_c)
        print("r_0 ", r_0)
        print("r_current", r_current)
        print("d", d)
        print("x_elem", x_elem)
        print("x_eq", x_eq)
        print("x_p", x_p)
        print("v_nucl", v_nucl)
        print("n_current", n_current)
        '''

        a = self.alpha_p(compound=compound, composition=current_solution, t=t)

        #print("a", a)

        try:
            part_1 = d / r_current
        except:
            part_1 = 0

        if part_1 < 0:
            part_1 = 0

        part_2 = (x_elem - (x_eq * math.exp(r_0 / r_current))) / ((a * x_p) - (x_eq * math.exp(r_0 / r_current)))

        r_pure = r_current + (part_1 * part_2 * d_tau) / 15.0

        r_corrected = (r_pure * n_current + d_n * r_0) / (n_current + d_n)

        return r_corrected

    def r_start(self, element="Nb", t=1350.0, current_solution={}):
        r_0 = self.r_0(element=element, t=t)
        x_elem = self.actual_x_concentration(element=element, current_solution=current_solution)
        x_eq = self.eq_x_concentration(element=element)
        return r_0 / (math.log((x_elem/x_eq), math.e)) + 1.0E-12

    def nucleation_rate_ostvald(self, compound="Nb_C_N", current_solution={}, t=1000.0, disl_dens=1.0, n_current=0.0,
                                r_current=1e-10):
        element = compound.split("_")[0]
        x_elem = self.actual_x_concentration(element=element, current_solution=current_solution)
        x_eq = self.eq_x_concentration(element=element)
        x_p = self.x_p_element(element=element)
        r_0 = self.r_0(element=element, t=t)
        d = self.diff_coef_eff(elemnt=element, r_current=r_current, t=t, disl_dens=disl_dens,
                               current_solution=current_solution)
        a = self.alpha_p(compound=compound, composition=current_solution, t=t)
        part_1 = 4.0 / 27.0
        part_2 = x_elem/((a*x_p) - x_eq)
        if r_current**3.0 != 0:
            part_3 = (r_0 * d)/(r_current**3.0)
        else:
            part_3 = 0
        part_4 = (r_0 * x_elem)/(r_current*(x_p - x_eq))
        part_5 = (3.0/(4.0*math.pi*(r_current**3.0))) / 1000.0 - n_current

        part_6 = 3.0*n_current
        '''
        print("---nucleation_rate_ostvald----")
        print("part_3", part_3)
        print("part_4", part_4)
        print("part_5", part_5)
        print("part_6", part_6)
        print("-----")
        '''


        return part_1 * part_2 * part_3 * ((part_4 * part_5) - part_6)

    def growth_rate_ostvald(self, element="Nb", t=1350.0, disl_dens=0.0, n_current=0.0,
                            r_current=1e-10, current_solution={}):
        d = self.diff_coef_eff(elemnt=element, r_current=r_current, t=t, disl_dens=disl_dens,
                               current_solution=current_solution)
        r_0 = self.r_0(element=element, t=t)
        x_elem = self.actual_x_concentration(element=element, current_solution=current_solution)
        x_eq = self.eq_x_concentration(element=element)
        x_p = self.x_p_element(element=element)
        compound = element + "_C_N"
        a = self.alpha_p(compound=compound, composition=current_solution, t=t)

        part_1 = 4.0 / 27.0
        part_2 = x_elem / ((a * x_p) - x_eq)
        part_3 = (r_0 * d) / (r_current ** 2.0)

        #if element == "Ti":
        #    print(r_0 * d, r_current ** 2.0)
        #    print(r_0, d, r_current)
        #    print(part_1, part_2, part_3)

        return part_1 * part_2 * part_3




class Solver:
    def __init__(self, parameters, composition, init_t=1350.0, D_0=83.0, e=0.35, v=2.0, base="Fe", lattice_type="FCC"):
        self.growth_handler = Growth(parameters, composition, init_t=init_t, base=base, lattice_type=lattice_type)

        self.rec_solver = rec_press.RecrystalizationSolver(parameters, composition)
        self.rec_solver.set_initial(D_0=D_0, e=e, v=v, T=init_t)

    @functools.lru_cache(maxsize=128, typed=False)
    def f_coars(self, r=1.0, r_c_eff=1.0):
        return 1.0 - math.erf(4.0 * ((r/r_c_eff) - 1))

    def check_complete(self, current_concentrations, eqiulibrium_concentrations, tolerance=0.01):
        max_diff = 0.0
        for i in current_concentrations:
            curr_diff = 2* (current_concentrations[i] - eqiulibrium_concentrations[i]) /\
                        (current_concentrations[i] + eqiulibrium_concentrations[i])
            if curr_diff > max_diff:
                max_diff = curr_diff
            if max_diff > tolerance:
                return False
        return True

    @functools.lru_cache(maxsize=128, typed=False)
    def calc_consumption(self, compound="Nb_C_N", r=0.0, n=0.0):
        #print("r= ", r)
        v = (4.0/3.0)*math.pi*n*(r**3.0)
        elements = compound.split("_")
        part_compound_1 = elements[0] + "_" + elements[1]
        part_compound_2 = elements[0] + "_" + elements[2]
        if self.growth_handler.parameters is not None:
            if "compounds_dens" in self.growth_handler.parameters:
                if (part_compound_1 in self.growth_handler.parameters["compounds_dens"] or
                    part_compound_2 in self.growth_handler.parameters["compounds_dens"]):
                    d_1 = 0.0
                    if part_compound_1 in self.growth_handler.parameters["compounds_dens"]:
                        d_1 = self.growth_handler.parameters["compounds_dens"][part_compound_1]
                    d_2 = 0.0
                    if part_compound_2 in self.growth_handler.parameters["compounds_dens"]:
                        d_2 = self.growth_handler.parameters["compounds_dens"][part_compound_2]
                    compound_composition = self.growth_handler.nucl_handler.compounds_compositions[compound]
                    d = d_1 * compound_composition[elements[1]] + d_2 * compound_composition[elements[2]]
                    m =  d * v
                    consumptions = {}
                    mr = 0
                    for i in elements:
                        mr += self.growth_handler.chem.atom_mass[i] * compound_composition[i]
                    if mr > 0:
                        for i in elements:
                            mass_part = (self.growth_handler.chem.atom_mass[i] * compound_composition[i]) / mr
                            consumptions[i] = mass_part * m
                        return consumptions
        return {}

    def calc_p_z(self, n_current={}, r_current={}, t=1000.0):
        g_gb = 1.41e6 * math.exp(-0.0117 * t)

        voluve_total = 0.0
        for j in n_current:
            #print("n_current", j, n_current[j])
            #print("r_current", j, r_current[j])
            voluve_total += (4.0 / 3.0) * math.pi * (r_current[j] ** 3.0) * n_current[j]

        r_sum = 0
        n_compounds = 0
        for j in r_current:
            r_sum += r_current[j]
            n_compounds += 1

        r_mean = 0
        if n_compounds > 0:
            r_mean = r_sum / float(n_compounds)

        p_z = 0.0

        #print('g_gb', g_gb)
        #print('voluve_total', voluve_total)
        #print('r_mean', r_mean)

        if r_mean > 0 and voluve_total > 0:
            p_z = (3.0 * g_gb * voluve_total) / r_mean

        return p_z

    def solve_isothermal(self, d_tau=0.01, max_x=0.999, max_step=100000, max_tau=1000, tolerance=0.00001):
        out_result = {}
        tau = 0.0
        t = self.growth_handler.current_t
        n_current = {}

        growth_rates = {}
        growth_rates_evalution = {}
        n_total = 0
        current_concentrations = self.growth_handler.chem.chemComposition
        equilibrium_concentrations = {}
        if "elements_in_solution" in self.growth_handler.equilibrium_data:
            equilibrium_concentrations = self.growth_handler.equilibrium_data["elements_in_solution"]
        else:
            return 0.0
        current_tau = 0.0
        current_step = 0

        r_current = {}

        #v_in_compounds = 0.0

        x_t = 0
        #x_t_prev = 0
        inner_m_p = 0

        c = self.growth_handler.chem.chemComposition["C"]
        #print(c)

        d_current = self.rec_solver.D_0

        p_z = self.calc_p_z(n_current=n_current, r_current=r_current, t=t)

        x_t, inner_m_p, disl_dens, d_current = self.rec_solver.x_t_calc(x_t, d_current, inner_m_p, c, current_tau,
                                                                        d_tau, p_z=p_z, n_total=n_total)
        n_init = self.growth_handler.nucl_handler.n_n(disl_dens=disl_dens, n=n_total)
        nucl_rates = self.growth_handler.nucl_handler.nucleation_rates(
            current_solution=current_concentrations, t=self.growth_handler.current_t, disl_dens=disl_dens, n=n_total)

        for i in nucl_rates:
            element = i.split("_")[0]
            t = self.growth_handler.current_t
            r_current[i] = self.growth_handler.r_start(element=element, t=t, current_solution=current_concentrations) #   .r_0(element=element, t=t)

        #print(current_concentrations)
        #print(equilibrium_concentrations)

        #print(self.check_complete(current_concentrations, equilibrium_concentrations, tolerance=tolerance))

        curent_elements_consumption = {}
        d_tau_initial = d_tau
        while  ((x_t < max_x) and (current_tau < max_tau) and (current_step < max_step)):
            # (not self.check_complete(current_concentrations, equilibrium_concentrations, tolerance=tolerance)
            n_init = self.growth_handler.nucl_handler.n_n(disl_dens=disl_dens, n=0)
            p_z = self.calc_p_z(n_current=n_current, r_current=r_current, t=t)
            #print(p_z)
            x_t, inner_m_p, disl_dens, d_current = self.rec_solver.x_t_calc(x_t, d_current, inner_m_p, c, current_tau,
                                                                            d_tau, p_z=p_z, n_total=n_total,
                                                                            n_init=n_init)

            nucl_rates = self.growth_handler.nucl_handler.nucleation_rates(
                current_solution=current_concentrations, t=self.growth_handler.current_t,
                disl_dens=disl_dens, n=n_total)

            for i in nucl_rates:
                if current_concentrations[i.split("_")[0]] > equilibrium_concentrations[i.split("_")[0]]:
                    if i in n_current:
                        t_s = self.growth_handler.nucl_handler.thermo_solver.find_solvus_carbonitride(i.split("_")[0])
                        if (current_concentrations[i.split("_")[0]] > equilibrium_concentrations[i.split("_")[0]]
                                and self.growth_handler.current_t < t_s):
                            n_current[i] += nucl_rates[i] * d_tau
                            n_total += n_current[i]
                    else:
                        n_current[i] = nucl_rates[i] * d_tau
                        n_total += n_current[i]

            #print("n_init", n_init)
            #print("n_total", n_total)

            #print(n_current)

            elements_consumptions = {}

            v_in_compounds = 0.0
            '''
            d_tau = self.estimate_d_tau(d_tau_initial=d_tau_initial, n_current=n_current,
                                        current_concentrations=current_concentrations,
                                        equilibrium_concentrations=equilibrium_concentrations,
                                        r_current=r_current, tolerance=tolerance, t=t, disl_dens=disl_dens,
                                        nucl_rates=nucl_rates, current_step=current_step)
            print('d_tau', d_tau)
            '''
            for i in n_current:
                element = i.split("_")[0]
                t_s = self.growth_handler.nucl_handler.thermo_solver.find_solvus_carbonitride(element=element)
                if (current_concentrations[element] > equilibrium_concentrations[element]
                    and self.growth_handler.current_t < t_s):
                    # and r_current[i] < self.growth_handler.r_c_eff(t=t, compound=i, current_solution=current_concentrations)):

                    gr_diff = 100.0
                    r_current_i = r_current[i]
                    r_current_i_2 = r_current[i]

                    r_prev = r_current_i
                    counter = 0

                    d_r = 0.0

                    if n_current[i] > 0:
                        '''
                        growth_rate_tmp = self.growth_handler.growth_rate(element=element, d_tau=d_tau, t=t, disl_dens=disl_dens,
                                                                          n_current=n_current[i], r_current=r_current_i,
                                                                          v_nucl=nucl_rates[i],
                                                                          current_solution=current_concentrations)
                        '''
                        r_current_i = self.growth_handler.new_radius(element=element, d_n=nucl_rates[i] * d_tau,
                                                                     d_tau=d_tau, t=t, disl_dens=disl_dens,
                                                                     n_current=n_current[i] - nucl_rates[i] * d_tau,
                                                                     r_current=r_current_i,
                                                                     v_nucl=nucl_rates[i],
                                                                     current_solution=current_concentrations)

                    else:
                        growth_rate_tmp = 0.0

                    #print("growth_rate_tmp_1", growth_rate_tmp)
                    # d_r = d_tau * growth_rate_tmp
                    '''
                    if growth_rate_tmp < 0:
                        d_r = 0.0
                        growth_rate_tmp = 0.0

                    if r_current_i <= 1.0001e-11 and i not in r_current:
                        d_r = d_tau * (growth_rate_tmp + growth_rate_prev) / 2.0
                        r_current_i = d_r
                    elif current_step < 1:
                        d_r = d_tau * (growth_rate_tmp + 0) / 2.0
                        r_current_i = r_current[i] + d_r
                    else:
                        d_r = d_tau * (growth_rate_tmp + growth_rate_prev) / 2.0
                        r_current_i = r_current[i] + d_r
                    '''

                    #growth_rate_prev = growth_rate_tmp

                    if n_current[i] > 0:
                        '''                        growth_rate_tmp = self.growth_handler.growth_rate(element=element, d_tau=d_tau, t=t, disl_dens=disl_dens,
                                                                          n_current=n_current[i], r_current=r_current_i,
                                                                          v_nucl=nucl_rates[i],
                                                                          current_solution=current_concentrations)
                        '''

                        r_current_i_2 = self.growth_handler.new_radius(element=element, d_n=nucl_rates[i] * d_tau,
                                                                     d_tau=d_tau, t=t, disl_dens=disl_dens,
                                                                     n_current=n_current[i] - nucl_rates[i] * d_tau,
                                                                     r_current=r_prev,
                                                                     v_nucl=nucl_rates[i],
                                                                     current_solution=current_concentrations)
                    '''
                    if growth_rate_tmp < 0:
                        d_r = 0.0
                        growth_rate_tmp = 0.0
                        growth_rates[i] = 0.0
                    '''

                    #print("growth_rate_tmp_2", growth_rate_tmp)

                    r_current_i = (r_current_i + r_current_i_2) / 2.0

                    growth_rates[i] = (r_current_i - r_current[i]) / d_tau

                    '''
                    if r_current_i <= 1.0001e-11 and i not in r_current:
                        #d_r = d_tau * (growth_rate_tmp + growth_rate_prev) / 2.0
                        growth_rates[i] = (r_current_i - r_current[i])/d_tau
                        r_current_i = d_r
                    elif current_step < 1:
                        d_r = d_tau * (growth_rate_tmp + 0) / 2.0
                        growth_rates[i] = (growth_rate_tmp + 0) / 2.0
                        r_current_i = r_current[i] + d_r
                    else:
                        d_r = d_tau * (growth_rate_tmp + growth_rate_prev) / 2.0
                        growth_rates[i] = (growth_rate_tmp + growth_rate_prev) / 2.0
                        r_current_i = r_current[i] + d_r
                    '''

                    #print("growth_rate_tmp_3", growth_rates[i])

                    if i in growth_rates_evalution:
                        growth_rates_evalution[i].append(growth_rates[i])
                    else:
                        growth_rates_evalution[i] = [growth_rates[i]]

                    if r_current_i > 0:
                        r_current[i] = r_current_i

                    #d_v_p_list[i] += (4.0 / 3.0) * math.pi * n_current[i] * r_current[i]

                    #print(n_current)
                    #print(r_current)

                    v_in_compounds += (4.0 / 3.0) * math.pi * n_current[i] * (r_current[i]**3.0)

                    cons_current = self.calc_consumption(compound=i, r=r_current[i], n=n_current[i])

                    for j in cons_current:
                        if j in elements_consumptions:
                            elements_consumptions[j] += cons_current[j]
                        else:
                            elements_consumptions[j] = cons_current[j]

                else:
                    # Ostwald riping
                    t = self.growth_handler.current_t
                    nucl_rates[i] = self.growth_handler.nucleation_rate_ostvald(compound=i, current_solution=current_concentrations,
                                                                                t=t, n_current=n_current[i], r_current=r_current[i])

                    # print("ostv", i, nucl_rates[i])
                    # if i == "Ti_C_N":
                    #     print(n_current[i] - nucl_rates[i] * d_tau)
                    #    print(nucl_rates[i])

                    if (n_current[i] - nucl_rates[i] * d_tau) > 0:

                        if nucl_rates[i] > 0:
                            n_current[i] -= nucl_rates[i] * d_tau
                            n_total -= nucl_rates[i] * d_tau
                        else:
                            n_current[i] += nucl_rates[i] * d_tau
                            n_total += nucl_rates[i] * d_tau

                        gr_diff = 100.0
                        r_current_i = r_current[i]
                        growth_rate_prev = 0.0
                        counter = 0

                        while gr_diff > 0.000001 and counter < 100:
                            counter += 1
                            growth_rate_tmp = self.growth_handler.growth_rate_ostvald(element=element, t=t, disl_dens=disl_dens,
                                                                n_current=n_current[i], r_current=r_current[i],
                                                                current_solution=current_concentrations)
                            #print("growth_rate_tmp_ostv", growth_rate_tmp)
                            # d_r = (d_r + d_tau * growth_rate_tmp) / 2.0
                            #if growth_rate_tmp != 0:
                            #    d_tau = r_current[i] / growth_rate_tmp / 10.0
                            d_r = d_tau * growth_rate_tmp

                            if r_current_i <= 1.0001e-10:
                                r_current_i = d_r
                            else:
                                r_current_i = r_current[i] + d_r

                            gr_diff = (2.0 * (growth_rate_prev - growth_rate_tmp) / (
                            growth_rate_prev + growth_rate_tmp)) ** 2.0

                            growth_rate_prev = growth_rate_tmp

                            # print("gr_diff ", gr_diff)

                        r_current[i] = r_current_i
                        growth_rates[i] = self.growth_handler.growth_rate_ostvald(element=element, t=t, disl_dens=disl_dens,
                                                                n_current=n_current[i], r_current=r_current[i],
                                                                current_solution=current_concentrations)

                        #print(i, current_tau, r_current[i])

                        if r_current[i] <= 1.0001e-10:
                            r_current[i] = d_tau * growth_rates[i]
                        else:
                            r_current[i] += d_tau * growth_rates[i]

            # print("t", current_tau)
            # print("r_NB", r_current["Nb_C_N"])
            # print("r_Ti", r_current["Ti_C_N"])

            if "x_t" not in out_result:
                out_result["x_t"] = {current_tau: x_t}
            else:
                out_result["x_t"][current_tau] = x_t
            if "growth_rates" not in out_result:
                out_result["growth_rates"] = {current_tau: copy.deepcopy(growth_rates)}
            else:
                out_result["growth_rates"][current_tau] = copy.deepcopy(growth_rates)
            if "n_current" not in out_result:
                out_result["n_current"] = {current_tau: copy.deepcopy(n_current)}
            else:
                out_result["n_current"][current_tau] = copy.deepcopy(n_current)
            if "r_current" not in out_result:
                out_result["r_current"] = {current_tau: copy.deepcopy(r_current)}
            else:
                out_result["r_current"][current_tau] = copy.deepcopy(r_current)
            if "current_concentrations" not in out_result:
                out_result["current_concentrations"] = {current_tau: copy.deepcopy(current_concentrations)}
            else:
                out_result["current_concentrations"][current_tau] = copy.deepcopy(current_concentrations)
            if "disl_dens" not in out_result:
                out_result["disl_dens"] = {current_tau: copy.deepcopy(disl_dens)}
            else:
                out_result["disl_dens"][current_tau] = copy.deepcopy(disl_dens)
            if "d_current" not in out_result:
                out_result["d_current"] = {current_tau: copy.deepcopy(d_current)}
            else:
                out_result["d_current"][current_tau] = copy.deepcopy(d_current)


            #print(x_t)
            #print(growth_rates["Ti_C_N"])
            #print(n_current["Ti_C_N"])
            #print(r_current["Ti_C_N"])
            #print(current_concentrations["Ti"])
            #print(self.growth_handler.r_c_eff(t=self.growth_handler.current_t, compound="Ti_C_N", current_solution=current_concentrations))
            #print(equilibrium_concentrations)
            #print(v_in_compounds)
            #print()

            aus_denf = self.growth_handler.solution_dens(t=self.growth_handler.current_t, composition=current_concentrations)

            masses_in_v = {}
            for j in current_concentrations:
                masses_in_v[j] = aus_denf * current_concentrations[j] / 100.0

            for j in masses_in_v:
                if j in elements_consumptions:
                    masses_in_v[j] -= elements_consumptions[j]
                    if masses_in_v[j] < 0:
                        masses_in_v[j] = equilibrium_concentrations[element] - (equilibrium_concentrations[element] / 1000.0)

            #print(masses_in_v)

            new_concentrations = {}
            mew_mass = 0.0
            for j in masses_in_v:
                mew_mass += masses_in_v[j]

            for j in masses_in_v:
                new_concentrations[j] = 100* masses_in_v[j] / mew_mass

            # print("new_concentrations", new_concentrations["Nb"])

            #print(new_concentrations)

            current_concentrations = new_concentrations

            new_disl_handler = rec_press.Dislocations(parameters=self.growth_handler.parameters, composition=current_concentrations)
            self.growth_handler.disl_handler = new_disl_handler
            self.growth_handler.nucl_handler.disl_handler = new_disl_handler
            self.rec_solver.dislocations_handler = new_disl_handler

            c = current_concentrations["C"]

            #print(self.f_coars(r=r_current["Ti_C_N"], r_c_eff=self.growth_handler.r_c_eff(t=self.growth_handler.current_t, compound="Ti_C_N", current_solution=current_concentrations)))
            #print(current_tau)

            current_tau += d_tau
            current_step += 1

        out_result["temperature"] = self.growth_handler.current_t
        return out_result

    @functools.lru_cache(maxsize=128, typed=False)
    def estimate_d_tau(self, d_tau_initial=0.01, n_current={}, current_concentrations={}, equilibrium_concentrations={},
                       r_current={}, tolerance=0.00001, t=1000, disl_dens=100000000, nucl_rates={}, current_step=0):
        current_d_tau = d_tau_initial
        d_tau = current_d_tau
        for i in n_current:
            element = i.split("_")[0]
            t_s = self.growth_handler.nucl_handler.thermo_solver.find_solvus_carbonitride(element=element)
            if (current_concentrations[element] > equilibrium_concentrations[element]
                    and self.growth_handler.current_t < t_s):

                gr_diff = 100.0
                r_current_i = r_current[i]
                growth_rate_prev = 0.0
                r_prev = 0.0
                counter = 0

                d_r = 0.0

                while counter < 2:
                    counter += 1
                    if n_current[i] > 0:
                        growth_rate_tmp = self.growth_handler.growth_rate(element=element, d_tau=d_tau, t=t, disl_dens=disl_dens,
                                                                          n_current=n_current[i], r_current=r_current_i,
                                                                          v_nucl=nucl_rates[i],
                                                                          current_solution=current_concentrations)
                    else:
                        growth_rate_tmp = 0.0
                    if growth_rate_tmp != 0:
                        d_tau = r_current_i / growth_rate_tmp / 2.0

                    if current_d_tau > d_tau:
                        current_d_tau = d_tau

                    #print('growth_rate_tmp', growth_rate_tmp)
                    #print('r_current_i', r_current_i)

                    d_r = d_tau * growth_rate_tmp

                    if r_current_i <= 1.0001e-10 and i not in r_current:
                        d_r = d_tau * (growth_rate_tmp + growth_rate_prev) / 2.0
                        r_current_i = d_r
                    elif current_step < 1:
                        d_r = d_tau * (growth_rate_tmp + growth_rate_prev) / 2.0
                        r_current_i = d_r
                    else:
                        r_current_i = r_current[i] + d_r

                    if growth_rate_prev + growth_rate_tmp != 0:
                        gr_diff = (2.0 * (growth_rate_prev - growth_rate_tmp) / (
                                    growth_rate_prev + growth_rate_tmp)) ** 2.0
                    else:
                        gr_diff = 0

                    # gr_diff = (2.0 * (r_prev - r_current_i) / (
                    #    r_prev + r_current_i)) ** 2.0

                    # print(growth_rate_tmp, r_current_i)

                    growth_rate_prev = growth_rate_tmp
                    r_prev = r_current_i
            else:
                # Ostwald riping
                t = self.growth_handler.current_t
                r_current_i = r_current[i]

                nucl_rates[i] = self.growth_handler.nucleation_rate_ostvald(compound=i,
                                                                            current_solution=current_concentrations,
                                                                            t=t, n_current=r_current_i,
                                                                            r_current=r_current[i])

                if (n_current[i] - nucl_rates[i] * d_tau) > 0:
                    n_current[i] -= nucl_rates[i] * d_tau

                    gr_diff = 100.0
                    r_current_i = r_current[i]
                    growth_rate_prev = 0.0
                    counter = 0

                    while counter < 2:
                        counter += 1
                        growth_rate_tmp = self.growth_handler.growth_rate_ostvald(element=element, t=t,
                                                                                  disl_dens=disl_dens,
                                                                                  n_current=n_current[i],
                                                                                  r_current=r_current[i],
                                                                                  current_solution=current_concentrations)

                        # d_r = (d_r + d_tau * growth_rate_tmp) / 2.0
                        if growth_rate_tmp != 0:
                            d_tau = r_current_i / growth_rate_tmp

                        if current_d_tau > d_tau:
                            current_d_tau = d_tau

                        d_r = d_tau * growth_rate_tmp

                        if r_current_i <= 1.0001e-10:
                            r_current_i = d_r
                        else:
                            r_current_i = r_current_i + d_r

                        gr_diff = (2.0 * (growth_rate_prev - growth_rate_tmp) / (
                                growth_rate_prev + growth_rate_tmp)) ** 2.0

                        growth_rate_prev = growth_rate_tmp

        return current_d_tau











'''
nucl_handler = Nucleation(parameters, chemComposition, init_t=1350.0)

a = nucl_handler.conpounds_molar_volumes()
print(a)

for i in nucl_handler.equilibrium_data:
    print(i, nucl_handler.equilibrium_data[i])

print(nucl_handler.calc_carbonitrides_stoichiometry())

print(nucl_handler.form_list_of_alpha_n())

print(nucl_handler.calc_td_forces(current_solution=nucl_handler.thermo_solver.chem_composition))

print(nucl_handler.g_c(t=1350.0, compound="Nb_C_N", current_solution=nucl_handler.thermo_solver.chem_composition))
print(nucl_handler.r_c(t=1350.0, compound="Nb_C_N", current_solution=nucl_handler.thermo_solver.chem_composition))
print(nucl_handler.n_c(t=1350.0, compound="Nb_C_N", current_solution=nucl_handler.thermo_solver.chem_composition))
print(nucl_handler.zeldovich_factor(t=1350.0, compound="Nb_C_N", current_solution=nucl_handler.thermo_solver.chem_composition))

print(nucl_handler.beta(t=1350.0, compound="Nb_C_N", current_solution=nucl_handler.thermo_solver.chem_composition))

print(nucl_handler.nucleation_rates(current_solution=nucl_handler.thermo_solver.chem_composition, t=1350.0, disl_dens=1e10, n=0))

print("growth_test")

growth_handler = Growth(parameters, chemComposition, init_t=1350.0)
print(growth_handler.solution_dens(t=1350.0))
print(growth_handler.austenite_molar_volume(t=1350.0))

print(growth_handler.eq_x_concentration())
print(growth_handler.actual_x_concentration(element="Nb", current_solution=growth_handler.chem.chemComposition))
print(growth_handler.x_p_element())
print(growth_handler.alpha_p(t=1350.0))

t_s = growth_handler.nucl_handler.thermo_solver.find_solvus_carbonitride(element="Nb")
print(t_s)
print(growth_handler.gamma_p_gamma(t_sovus=t_s))
print(growth_handler.r_0(t=1350.0))
print(growth_handler.diff_coef(t=1350.0))
print(growth_handler.disl_diff_coef(t=1350.0))
print(growth_handler.diff_coef_eff(t=1350.0, disl_dens=1892308866164430.5, current_solution=growth_handler.chem.chemComposition))

solve_handler = Solver(parameters, chemComposition, init_t=1350.0)

print(solve_handler.f_coars(r=10.0, r_c_eff=1.0))
print(solve_handler.solve_isothermal(max_step=2500))
'''

#self.chem.chemComposition





