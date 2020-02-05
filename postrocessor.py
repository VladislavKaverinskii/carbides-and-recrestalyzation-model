# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import json
import data_loader
import xlwt


class Postprocessor:

   def __init__(self, path="terms_description.json"):
       self.dataset = None
       self.initial = None
       self.possible_to_plot = []
       self.possible_to_diagram = []
       self.prepeared_data = None
       self.terms_description = data_loader.load_terms_description(path=path)

   def clear_data(self):
       self.possible_to_plot = []
       self.possible_to_diagram = []
       self.prepeared_data = None

   def load_dataset(self, file):
       self.clear_data()
       f = open(file, 'r')
       raw_data = json.loads(f.read())
       if "result" in raw_data:
           self.dataset = raw_data["result"]
       if "initial" in raw_data:
           self.initial = raw_data["initial"]
       f.close()


   def prepear_data(self):
       output = {}
       self.clear_data()
       if self.dataset is not None:
           for key in self.dataset:
               output[key] = {"name": key}
               if isinstance(self.dataset[key], dict):
                   try:
                       independent = []
                       depandents = {}
                       for key_2 in self.dataset[key]:
                           try:
                               independent.append(float(key_2))
                               if isinstance(self.dataset[key][key_2], dict):
                                   for key_3 in self.dataset[key][key_2]:
                                       if key_3 in depandents:
                                           depandents[key_3].append(self.dataset[key][key_2][key_3])
                                       else:
                                           depandents[key_3] = [self.dataset[key][key_2][key_3]]
                               elif isinstance(self.dataset[key][key_2], float) or isinstance(self.dataset[key][key_2], int):
                                   if self.terms_description[key]["y"] in depandents:
                                       depandents[self.terms_description[key]["y"]].append(self.dataset[key][key_2])
                                   else:
                                       depandents[self.terms_description[key]["y"]] = [self.dataset[key][key_2]]
                               if key not in self.possible_to_plot:
                                   self.possible_to_plot.append(key)
                           except ValueError:
                               if key not in self.possible_to_diagram:
                                   self.possible_to_diagram.append(key)
                               independent.append(key_2)
                               if "y" in depandents:
                                   depandents["y"].append(self.dataset[key][key_2])
                               else:
                                   depandents["y"] = [self.dataset[key][key_2]]

                       output[key]["x"] = independent
                       output[key]["y"] = depandents

                   except:
                       pass
               else:
                   output[key]["data"] = self.dataset[key]
       self.prepeared_data = output
       return output

   def make_plot(self, parameter_key):
       if parameter_key in self.possible_to_plot and self.prepeared_data is not None:
           data = self.prepeared_data[parameter_key]
           plot_label = self.terms_description[parameter_key]["name"]
           x = data["x"]
           x_label = self.terms_description[parameter_key]["x"]
           y_label = self.terms_description[parameter_key]["y"]
           for data_list in data["y"]:
               if data_list == "y":
                   line_name = self.terms_description[parameter_key]["y"]
               else:
                   line_name = data_list
               if line_name.lower() != "fe":
                   y = data["y"][data_list]
                   plt.plot(x, y, label=line_name)

           plt.xlabel(x_label)
           plt.ylabel(y_label)
           plt.title(plot_label)
           plt.grid(True)
           plt.rc('grid', linestyle="-", color='black')
           plt.legend()
           plt.show()
       if parameter_key in self.possible_to_diagram and self.prepeared_data is not None:
           data = self.prepeared_data[parameter_key]
           plot_label = self.terms_description[parameter_key]["name"]
           x = data["x"]
           counter = 0
           deleted_key = False
           while counter < len(x):
               if x[counter].lower() == "fe":
                   del x[counter]
                   deleted_key = True
                   break
               counter += 1

           y_label = self.terms_description[parameter_key]["y"]
           x_label = self.terms_description[parameter_key]["x"]
           for data_list in data["y"]:
               if data_list == "y":
                   line_name = self.terms_description[parameter_key]["y"]
               else:
                   line_name = data_list
               if line_name.lower() != "fe":
                   y = data["y"][data_list]
                   if deleted_key:
                       del y[counter]
                   plt.bar(x, y, width=0.2, alpha=0.7, label=line_name, zorder=2)

           plt.xlabel(x_label)
           plt.ylabel(y_label)
           plt.title(plot_label)
           plt.legend()
           plt.show()

   def save_to_xls(self, fail_name="out.xls", every=1):
       existing_sheets = []
       if self.prepeared_data is None:
           self.prepear_data()
       book = xlwt.Workbook(encoding="utf-8")
       for parameter_key in self.prepeared_data:
           if len(self.terms_description[parameter_key]["name"]) < 10:
               if self.terms_description[parameter_key]["name"] not in existing_sheets:
                   new_sheet = book.add_sheet(self.terms_description[parameter_key]["name"])
                   existing_sheets.append(self.terms_description[parameter_key]["name"])
               else:
                   counter = 0
                   new_sheet_name = self.terms_description[parameter_key]["name"] + str(counter)
                   while new_sheet_name in existing_sheets:
                       new_sheet_name = self.terms_description[parameter_key]["name"] + str(counter)
                       counter += 1
                   new_sheet = book.add_sheet(new_sheet_name)
                   existing_sheets.append(new_sheet_name)
           else:
               if self.terms_description[parameter_key]["name"][:10] not in existing_sheets:
                   new_sheet = book.add_sheet(self.terms_description[parameter_key]["name"][:10])
                   existing_sheets.append(self.terms_description[parameter_key]["name"][:10])
               else:
                   counter = 0
                   new_sheet_name = self.terms_description[parameter_key]["name"][:10] + str(counter)
                   while new_sheet_name in existing_sheets:
                       new_sheet_name = self.terms_description[parameter_key]["name"][:10] + str(counter)
                       counter += 1
                   new_sheet = book.add_sheet(new_sheet_name)
                   existing_sheets.append(new_sheet_name)
           new_sheet.write(0, 0, self.terms_description[parameter_key]["x"])
           if "y" in self.prepeared_data[parameter_key]:
               counter = 1
               for sub_parameter in self.prepeared_data[parameter_key]["y"]:
                   if sub_parameter == "y":
                       new_sheet.write(0, counter, self.terms_description[parameter_key]["y"])
                   else:
                       new_sheet.write(0, counter, str(sub_parameter))
                   counter += 1
               counter = 1
               for value in self.prepeared_data[parameter_key]["x"]:
                   if (counter-1)%every == 0:
                       new_sheet.write(int((counter-1)/every + 1), 0, value)
                   counter += 1
               counter = 1
               for sub_parameter in self.prepeared_data[parameter_key]["y"]:
                   counter_2 = 1
                   for value in self.prepeared_data[parameter_key]["y"][sub_parameter]:
                       if (counter_2 - 1) % every == 0:
                           new_sheet.write(int((counter_2-1)/every + 1), counter, value)
                       counter_2 += 1
                   counter += 1
           elif "data" in self.prepeared_data[parameter_key]:
               new_sheet.write(1, 1, self.terms_description[parameter_key]["name"])
               new_sheet.write(1, 2, self.prepeared_data[parameter_key]["data"])
       book.save(str(fail_name))








my_Postprocessor = Postprocessor()

my_Postprocessor.load_dataset("carbides_kinetics_result.json")
print(my_Postprocessor.dataset)
print(my_Postprocessor.prepear_data())

my_Postprocessor.make_plot("x_t")
my_Postprocessor.make_plot("growth_rates")
my_Postprocessor.make_plot("n_current")
my_Postprocessor.make_plot("r_current")
my_Postprocessor.make_plot("current_concentrations")
my_Postprocessor.make_plot("disl_dens")
my_Postprocessor.make_plot("d_current")

my_Postprocessor.save_to_xls()

print()

my_Postprocessor.load_dataset("recrystalization_result.json")
d = my_Postprocessor.prepear_data()
for i in d:
    print(i)
    print(d[i])

print()

my_Postprocessor.make_plot("grain_d")
#my_Postprocessor.save_to_xls()

my_Postprocessor.load_dataset("thermodynamic_result.json")
d = my_Postprocessor.prepear_data()
for i in d:
    print(i)
    print(d[i])

my_Postprocessor.make_plot('elements_consumption')
my_Postprocessor.make_plot('percents_of_compounds')
my_Postprocessor.make_plot('elements_in_solution')

# my_Postprocessor.save_to_xls()









