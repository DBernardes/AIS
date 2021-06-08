# -*- coding: utf-8 -*-
"""Flux Calculation class tests.

This script tests the operation of the Background Image Class.

Created on Thu Apr 22 13:44:35 2021

@author: denis
"""


from code.BGI.BGI import Background_Image
import pytest

dic = {"em_mode": 0, "em_gain": 1, "binn": 1, "t_exp": 1, "preamp": 1, "hss": 1}


@pytest.fixture
def bgi():
    return Background_Image(
        ccd_operation_mode=dic,
        ccd_gain=3,
        dark_current=1e-5,
        read_noise=6.3,
        bias_level=500,
    )


# ------------------------ Initialize the class --------------------------


def test_em_gain(bgi):
    assert bgi.em_gain == 1


def test_bin(bgi):
    assert bgi.binn == 1


def test_t_exp(bgi):
    assert bgi.t_exp == 1


def test_ccd_gain(bgi):
    assert bgi.ccd_gain == 3


def test_bias_level(bgi):
    assert bgi.bias_level == 500


# ----------------------- Calculate Background Image -------------------------


def test_create_background_image(bgi):
    bgi.create_background_image(100)


def test_create_bias_image(bgi):
    bgi.create_bias_image()


def test_create_dark_image(bgi):
    bgi.create_dark_image()
