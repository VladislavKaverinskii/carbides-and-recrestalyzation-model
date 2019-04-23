import data_loader
import rec_press
import carbides
import thermodytamic
import json


class MainSolver:

    def __init__(self):
        self.task = None
        self.parameters = None
        self.thermodynamic_result = None
        self.carbides_result = None
        self.recrystalization_result = None


    def read_data(self, parameters = "parameters.json", task = "task.json"):
        self.task = data_loader.load_task(path=task)
        self.parameters = data_loader.load_parameters(path=parameters)
        return (self.task, self.parameters)


    def calculate_carbonitrides_equilibrium(self):
        if self.task is None and self.parameters is None:
            self.task, self.parameters = self.read_data(parameters = "parameters.json", task = "task.json")
        if self.task is not None:
            if self.task.base == "Fe":
                thermo_solver = thermodytamic.CarbonitrideThermodynamicsInSteel(chem_composition=self.task.chemComposition,
                                                                                td_params=self.parameters["td_params"],
                                                                                t=self.task.init_t)
                if thermo_solver is not None:
                    self.thermodynamic_result = thermo_solver.solve()
                    return self.thermodynamic_result
        return None


    def carbonitride_precipitation_simulation(self):
        if self.task is None and self.parameters is None:
            self.task, self.parameters = self.read_data(parameters = "parameters.json", task = "task.json")
        if self.task is not None:
            if self.task.base == "Fe":
                solve_handler = carbides.Solver(self.parameters, self.task.chemComposition, init_t=self.task.init_t,
                                                D_0=self.task.D_0, e=self.task.e_deform, v=self.task.v_deform,
                                                base=self.task.base, lattice_type=self.task.lattice_type)
                if solve_handler is not None:
                    self.carbides_result = solve_handler.solve_isothermal(d_tau=self.task.d_tau, max_x=self.task.max_x, max_step=self.task.max_step,
                                                          max_tau=self.task.max_tau, tolerance=self.task.tolerance)
                    return self.carbides_result
        return None


    def recrystalization_simulation(self):
        if self.task is None and self.parameters is None:
            self.task, self.parameters = self.read_data(parameters = "parameters.json", task = "task.json")
        if self.task is not None:
            if self.task.base == "Fe":
                solver = rec_press.RecrystalizationSolver(self.parameters, self.task.chemComposition)
                solver.set_initial(D_0=self.task.D_0, e=self.task.e_deform, v=self.task.v_deform, T=self.task.init_t)
                self.recrystalization_result = solver.solve_austenite(n_step_max=self.task.max_step, d_tau=self.task.d_tau,
                                                                      max_tau=self.task.max_tau, max_x=self.task.max_x)
                return self.recrystalization_result
        return None


    def save_thermodynamic_result(self, out_file="thermodynamic_result.json"):
        if self.thermodynamic_result is None:
            self.thermodynamic_result = self.calculate_carbonitrides_equilibrium()
        f = open(out_file, 'w')
        f.write(json.dumps({"initial": self.task.serialize(), "result": self.thermodynamic_result}))
        f.close()


    def save_carbides_result(self, out_file="carbides_kinetics_result.json"):
        if self.carbides_result is None:
            self.carbides_result = self.carbonitride_precipitation_simulation()
        f = open(out_file, 'w')
        f.write(json.dumps({"initial": self.task.serialize(), "result": self.carbides_result}))
        f.close()


    def save_recrystalization_result(self, out_file="recrystalization_result.json"):
        if self.recrystalization_result is None:
            self.recrystalization_result = self.carbonitride_precipitation_simulation()
        f = open(out_file, 'w')
        f.write(json.dumps({"initial": self.task.serialize(), "result": self.recrystalization_result}))
        f.close()


if __name__ == "__main__":
    solver = MainSolver()
    solver.read_data()
    print(solver.calculate_carbonitrides_equilibrium())
    solver.save_thermodynamic_result()
    print(solver.carbonitride_precipitation_simulation())
    solver.save_carbides_result()
    print(solver.recrystalization_simulation())
    solver.save_recrystalization_result()








