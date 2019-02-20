# -*- coding: utf-8 -*-

import math

R = 8.3144598

A_EXP_COEF = 0.00007076

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
    "C": 0.1,
    "Mn": 0.65,
    "Si": 0.2,
    "Mo": 0.02,
    "Ti": 0.006,
    "V": 0.005,
    "Nb": 0.005
}

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

parameters = {
    "Qd_params": Qd_params,
    "A_params": A_params,
    "pic_def_params_n": pic_def_params_n,
    "pic_def_params_A": pic_def_params_A
}


class EmpiricParams():
    def B(self, Z=1.0, A=1.0):
        if (Z!=0) and (A!=0):
            tmp = 9.5326+0.6196*math.log((Z/A), math.e)
            return tmp*tmp*1E6
        return 0

    def C(self, Z=1.0, A=1.0):
        if A != 0:
            return 3.9202*(math.pow((Z/A), 0.0592))
        return 0

    def m(self, Z=1.0, A=1.0):
        if A != 0:
            if (Z/A) >=0:
                return 0.3449*math.exp(0.0139*math.sqrt(Z/A))
        return 0

    def B_(self, Z=1.0, A=1.0):
        if A != 0:
            return 26.0310*math.pow((Z/A), 0.1351)*1E6
        return 0

    def k(self, Z=1.0, A=1.0):
        if Z != 0:
            if (Z/A) > 0:
                return 0.5974*math.exp((1.2333*math.sqrt(A/Z)))
        return 0

    def m_(self, Z=1.0, A=1.0):
        if A != 0:
            if (Z/A) >=0:
                return 1.0901*math.exp(0.264*math.sqrt(Z/A))
        return 0

    def sum_param(self, params, chem_composition):
        n = 0
        if "intercept" in params:
            n = float(params["intercept"])
        for i in chem_composition:
            if i in params:
                n += float(params[i])*float(chem_composition[i])
        return n


class FlowStress():

    def __init__(self, parameters={}, chem_composition={}):
        self.qd_params = parameters['Qd_params']
        self.chem_composition = chem_composition
        self.A_params = parameters['A_params']
        self.pic_def_params_n = parameters['pic_def_params_n']
        self.pic_def_params_A = parameters['pic_def_params_A']

        self.empiric = EmpiricParams()


    def __get_se__(self, e=0, B=0, C=0, m=0):
        return B*pow((1-math.exp(-C*e)), m)


    def get_Qd(self):
        sum = 0
        if "intercept" in self.qd_params:
            sum = self.qd_params["intercept"]
        for i in self.chem_composition:
            coef = 0
            if (i + "_coef") in self.qd_params:
                try:
                    coef = float(self.qd_params[i+"_coef"])
                except:
                    coef = 0
            if coef != 0:
                index = 1
                if (i+"_index") in self.qd_params:
                    try:
                        index = float(self.qd_params[i+"_index"])
                    except:
                        index = 1
                if index != 1:
                    sum += coef*pow(self.chem_composition[i], index)
                else:
                    sum += coef * self.chem_composition[i]
        return sum

    def get_Z(self, v=0.0, T=1.0):
        if v != 0:
            try:
                return v*math.exp(self.get_Qd()/ (R*T))
            except:
                return 0
        else:
            return 0

    def get_A(self):
        factor = 0
        if "intercept" in self.A_params:
            factor = float(self.A_params["intercept"])
        for i in self.chem_composition:
            if i in self.A_params:
                factor += self.A_params[i]*self.chem_composition[i]
        if factor != 0:
            return factor*math.exp(A_EXP_COEF * self.get_Qd())
        else:
            return 0

    def e_p(self, D_0=0.0, v=0.0, T=1.0):
        e_0 = self.empiric.sum_param(self.pic_def_params_n, self.chem_composition)
        d = math.pow(D_0, 0.195)
        Z_A = self.get_Z(v=v, T=T)/self.get_A()
        n_2 = self.empiric.sum_param(self.pic_def_params_A, self.chem_composition)
        return  (e_0 * d * pow(Z_A , n_2))

    def f_DRX(self, D_0=0.0, e=0.0, v=1.0, T=1.0):
        e_p = self.e_p(D_0=D_0, v=v, T=T)
        k = self.empiric.k(Z=self.get_Z(v=v, T=T), A=self.get_A())
        m = self.empiric.m_(Z=self.get_Z(v=v, T=T), A=self.get_A())
        if e > 0.95*e_p:
            a1 = (e-(0.95*e_p))/e_p
            a2 =  pow(a1, m).real
            return (1 - math.exp(-k*a2))
        else:
            return 0

    def rec_sigma(self, D_0=0.0, e=0.0, v=1.0, T=1.0):
        b = self.empiric.B_(Z=self.get_Z(v=v, T=T), A=self.get_A())
        f = self.f_DRX(D_0=D_0, e=e, v=v, T=T)
        return b * f

    def ret_sigma(self, e=0.0, v=1.0, T=1.0):
        z = self.get_Z(v=v, T=T)
        a = self.get_A()
        c = self.empiric.C(Z=z, A=a)
        m = self.empiric.m(Z=z, A=a)
        b = self.empiric.B(Z=z, A=a)
        return b * pow((1 - math.exp(-c * e)), m)

    def calcFlowStress(self, D_0=0.0, e=0.0, v=1.0, T=1.0):
        return self.ret_sigma(e=e, v=v, T=T) - self.rec_sigma(D_0=D_0, e=e, v=v, T=T)

    def d_DRX(self, v=1.0, T=1.0):
        z = self.get_Z(v=v, T=T)
        a = self.get_A()
        return 38.26*pow((z/a), -0.08)

'''
fs = FlowStress(parameters=parameters, chem_composition=chemComposition)


s = fs.rec_sigma(D_0=50, e=0.5, v=0.00011, T=1200.0)
print(s)
s_2 = fs.ret_sigma(e=0.5, v=0.00011, T=1200.0)
print(s_2)
f_s = fs.calcFlowStress(D_0=50, e=0.5, v=0.00011, T=1200.0)
print(f_s)
drx = fs.d_DRX(v=0.00011, T=1200.0)
print(drx)
'''



























