
import numpy as np
import physicalconstants as const



class Gas(object):

    def __init__(self, rho, lmbda, mu, Cp, beta):

        self.rho = rho  # kg/m3
        self.lmbda = lmbda  # W/m/K
        self.mu = mu # kg/m/s
        self.Cp = Cp # J/kg/K
        self.beta = beta # 1/K

    @property
    def Prandtl(self):
        Pr = self.Cp * self.mu / self.lmbda
        return Pr




class Gap(object):
    """ Gaz gap layer object
        Compute the thermal conductance due to convection
        as a function of deltaT
    """

    def __init__(self, gas, thickness):
        self.gas = gas
        self.w = thickness

    def Grashoff(self, deltaT):
        deltaT = np.abs( deltaT )
        gas = self.gas
        Gr = const.g * gas.beta*gas.rho**2*self.w**3 / gas.mu**2 * deltaT
        return Gr

    def Nusselt(self, deltaT):
        """ Nusselt number
            Accept float or nd-array
        """
        deltaT = np.asarray( deltaT )

        GrPr = self.gas.Prandtl * self.Grashoff( deltaT )

        condlist = [ GrPr < 5e3,
                     np.logical_and( GrPr>5e3 , GrPr<6e4 ),
                     np.logical_and( GrPr>6e6 , GrPr<1.5e5 ),
                     GrPr > 1.5e5  ]
        funclist = [ 1,
                     lambda x: 0.0429 * x**0.37 ,
                     lambda x: 0.43 * x**0.16 ,
                     lambda x: 0.0354 * x**0.37  ]

        Nu = np.piecewise(GrPr, condlist, funclist)

        return Nu

    def h(self, deltaT):
        """ conductance in W/m2/K
        """
        h = self.gas.lmbda / self.w * self.Nusselt( deltaT )
        return h



class MaterialLibrary(object):
    """ common physical properties for materials
    """

    @property
    def air(self):
        return Gas( 1.29, 2.5e-2, 1.86e-5, 1.005e3, 3.67e-3 )

matlib = MaterialLibrary()
