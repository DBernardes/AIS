# -*- coding: utf-8 -*-

"""Channel Creator Package.

This is the package of the channel creators. The AIS is based on the Factory
Method to implement the correct star flux calculation as a function of the
optical path of each channel. This package has as abstract class that will
be provided to the Point_Spread_Function Class. Then, the Point_Spread_Function
class will call the correct concrete channel creator for the respective desired
image.
"""


class Abstract_Channel_Creator:

    def __init__(self):
        pass

    def factory_method(self):
        pass

    def calc_star_flux(self):
        pass


class Concrete_Channel_1(Abstract_Channel_Creator):

    def factory_method(self):
        pass


class Concrete_Channel_2(Abstract_Channel_Creator):

    def factory_method(self):
        pass


class Concrete_Channel_3(Abstract_Channel_Creator):

    def factory_method(self):
        pass


class Concrete_Channel_4(Abstract_Channel_Creator):

    def factory_method(self):
        pass
