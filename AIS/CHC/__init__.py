# -*- coding: utf-8 -*-

"""Channel Creator Package.

This package creates an abstract channel as well as the four concrete channels
that simulates the 4 SPARC4 channels. This abstract channel is provided to the
Point_Spread_Function Class, that calls the correct concrete channel based on
the desired image.
"""
from .channel_creator import (Abstract_Channel_Creator,
                              Concrete_Channel_1,
                              Concrete_Channel_2,
                              Concrete_Channel_3,
                              Concrete_Channel_4)
