# -*- coding: utf-8 -*-

import math
import decimal
import functools

R = 8.3144598

A_EXP_COEF = 0.00007076


class EmpiricParams():
    @functools.lru_cache(maxsize=128, typed=False)
    def B(self, Z=1.0, A=1.0):
        if (Z!=0) and (A!=0) and (Z/A > 0):
            tmp = 9.5326+0.6196*math.log((Z/A), math.e)
            return tmp*tmp*1E6
        return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def C(self, Z=1.0, A=1.0):
        if A != 0 and Z/A >= 0:
            return 3.9202*(math.pow((Z/A), 0.0592))
        return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def m(self, Z=1.0, A=1.0):
        # print("A", A)
        # print("Z", Z)
        if A != 0:
            if (Z/A) >= 0:
                try:
                    return 0.3449*math.exp(0.0139*math.sqrt(Z/A))
                except OverflowError:
                    return 1
        return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def B_(self, Z=1.0, A=1.0):
        if A != 0 and Z/A > 0:
            return 26.0310*math.pow((Z/A), 0.1351)*1E6
        return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def k(self, Z=1.0, A=1.0):
        if Z != 0:
            if (Z/A) > 0:
                try:
                    return 0.5974*math.exp((1.2333*math.sqrt(A/Z)))
                except OverflowError:
                    return 1
        return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def m_(self, Z=1.0, A=1.0):
        if A != 0:
            if (Z/A) >=0:
                try:
                    return 1.0901*math.exp(0.264*math.sqrt(Z/A))
                except OverflowError:
                    return 1
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

    @functools.lru_cache(maxsize=128, typed=False)
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

    @functools.lru_cache(maxsize=128, typed=False)
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
        if Z_A > 0:
            return  (e_0 * d * pow(Z_A , n_2))
        else:
            return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def f_DRX(self, D_0=0.0, e=0.0, v=1.0, T=1.0):
        e_p = self.e_p(D_0=D_0, v=v, T=T)
        k = self.empiric.k(Z=self.get_Z(v=v, T=T), A=self.get_A())
        m = self.empiric.m_(Z=self.get_Z(v=v, T=T), A=self.get_A())
        #print('e', e)
        #print('e_p', e_p)
        if e > 0.95*e_p and e_p != 0:
            a1 = (e-(0.95*e_p))/e_p
            a2 = pow(a1, m).real
            return (1 - math.exp(-k*a2))
        else:
            return 0

    @functools.lru_cache(maxsize=128, typed=False)
    def rec_sigma(self, D_0=0.0, e=0.0, v=1.0, T=1.0):
        b = self.empiric.B_(Z=self.get_Z(v=v, T=T), A=self.get_A())
        f = self.f_DRX(D_0=D_0, e=e, v=v, T=T)
        return b * f

    @functools.lru_cache(maxsize=128, typed=False)
    def ret_sigma(self, e=0.0, v=1.0, T=1.0):
        z = self.get_Z(v=v, T=T)
        a = self.get_A()
        c = self.empiric.C(Z=z, A=a)
        m = self.empiric.m(Z=z, A=a)
        b = self.empiric.B(Z=z, A=a)
        return b * pow((1 - math.exp(-c * e)), m)

    @functools.lru_cache(maxsize=128, typed=False)
    def calcFlowStress(self, D_0=0.0, e=0.0, v=1.0, T=1.0):
        return self.ret_sigma(e=e, v=v, T=T) - self.rec_sigma(D_0=D_0, e=e, v=v, T=T)

    @functools.lru_cache(maxsize=128, typed=False)
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



























