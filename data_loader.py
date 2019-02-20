# -*- coding: utf-8 -*-

import json

td_params = {
    "eq_constant_params": None,
    "vagner_interactions": None,
    "energetic_ponomarenko_params": None,
    "solute_energies": None,
    "formation_energies": None,
    "melting_temperatures": None,
    "h_melting": None
}


parameters = {
    "Qd_params": None,
    "A_params": None,
    "pic_def_params_n": None,
    "pic_def_params_A": None,
    "Qsd_parameters": None,
    "return_parametrs": None,
    "diffusion_parameters": None,
    "dislocations_diffusion": None,
    "compounds_dens": None,
    "td_params": td_params,
    "elements_dens": None,
    "surface_energy": None
}


class Task:
    def __init__(self):
        self.base = ""
        self.init_t = 0.0
        self.D_0 = 0.0
        self.e_deform = 0.0
        self.v_deform = 0.0
        self.chemComposition = {}
        self.lattice_type = ""
        self.max_step = 0
        self.r_init = 10.0
        self.d_tau = 0.01
        self.max_x = 0.999
        self.max_tau = 0
        self.tolerance = 0.00001


def load_parameters(path="parameters.json"):
    raw_json = json.loads(open(path).read())

    td_params = {
        "eq_constant_params": raw_json["eq_constant_params"],
        "vagner_interactions": raw_json["vagner_interactions"],
        "energetic_ponomarenko_params": raw_json["energetic_ponomarenko_params"],
        "solute_energies": raw_json["solute_energies"],
        "formation_energies": raw_json["formation_energies"],
        "melting_temperatures": raw_json["melting_temperatures"],
        "h_melting": raw_json["h_melting"]
    }

    parameters = {
        "Qd_params": raw_json["Qd_params"],
        "pic_def_params_n": raw_json["pic_def_params_n"],
        "pic_def_params_A": raw_json["pic_def_params_A"],
        "Qsd_parameters": raw_json["Qsd_parameters"],
        "return_parametrs": raw_json["return_parametrs"],
        "diffusion_parameters": raw_json["diffusion_parameters"],
        "dislocations_diffusion": raw_json["dislocations_diffusion"],
        "compounds_dens": raw_json["compounds_dens"],
        "td_params": td_params ,
        "elements_dens": raw_json["elements_dens"],
        "surface_energy": raw_json["surface_energy"]
    }

    return parameters


def load_task(path="task.json"):
    raw_json = json.loads(open(path).read())
    current_task = Task()

    current_task.base = raw_json["base"]
    current_task.init_t = raw_json["init_t"]
    current_task.D_0 = raw_json["D_0"]
    current_task.e_deform = raw_json["e_deform"]
    current_task.v_deform = raw_json["v_deform"]
    current_task.chemComposition = raw_json["chemComposition"]
    current_task.lattice_type = raw_json["lattice_type"]
    current_task.max_step = raw_json["max_step"]
    current_task.r_init = raw_json["r_init"]
    current_task.d_tau = raw_json["d_tau"]
    current_task.max_x = raw_json["max_x"]
    current_task.max_tau = raw_json["max_tau"]
    current_task.tolerance = raw_json["tolerance"]

    return current_task






