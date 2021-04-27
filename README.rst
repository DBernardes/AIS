|Unit Tests| |Documentation|

Introduction
============

This repo contains the software of the Artificial Images Simulator (AIS). The AIS was developed using python 3 language to 
create cubes of artificial star images, simulating those images that would be acquired by using the SPARC4 [#SPARC4]_ CCD cameras 
in astronomical observations. To create the images, the AIS models the star flux distribution as a 2D-Gaussian Distribution. 
This result is added to a noise image, created based on the noise information of the SPARC4 CCDs, as a function of their operation mode. 
Figure below presents an example of an image created by the AIS. This page explains the step-by-step procedure of the AIS to create the images. 
Also, it is presented a simple execution example. 

.. image:: /docs/images/artificial_star.png   
   :alt: Artificial star image
   :align: center


Software Description
====================

This section presents the procedure used the AIS to generate the cubes of artifical images. Initially, its operation can be divided into
the parameters configuration step, the creation of a noise image step and, the calculation of the star flux step. In the parameters configuration step, 
the AIS will obtain the gain, the dark current noise, and the read noise of the CCD, as a function of its operation mode. The CCD gain is the 
conversion factor between the acquired photoelectrons and the Analogical-to-Ditial Unit (ADU), and it is obtained through the camera datasheet. 
The dark current noise is the number of thermoelectrons per pixel created by the dark current of the CCD as a function of its respective temperature. 
The read noise is a fluctuation in the value measured in each pixel of the CCD resulting from the image readout process. The calculation of the dark 
current noise and the read noise of the camera are based on the characterization of the SPARC4 CCDs, presented by [#Bernardes_2018]_.

Then, it is created the noise image. To accomplish this, initially, it is calculated the background level (:math:`B_L`) of the image in ADU, as follows

.. image:: /docs/images/back_ground_level.png   
   :alt: Noise level equation
   :align: center
	

where B represents the bias level of the CCD in ADU; the :math:`S_{dc}` is the mean of the thermoelectrons created by the CCD dark current, 
in e-/pix; the :math:`S_{sky}` is the sky flux, in photons/pix; the :math:`G_{em}` is the Electron Multiplying Gain of the CCD; the B<sub>in</sub> 
is the number of pixels binned in the reading process; the G is the gain of the CCD. Once obtained the background level, it is calculated the image 
noise (N), in ADU. This noise is composed by the contributions of the sky noise, the dark current noise, and the read noise as follows

.. image:: /docs/images/noise_level.png   
   :alt: Noise level equation
   :align: center

where :math:`N_F` is an extra noise factor of the use of the Electron Multiplying mode, and it equals 1.4. For the Conventional Mode, 
:math:`N_F` = 1; :math:`\sigma_{ADU}` is the read noise in ADU. With these results, it is created the noise image, with pixels values given 
by a gaussian distribution centered in :math:`B_L` and a standard deviation of N.

The next step is the calculation of the star flux distribution. To accomplish this, it is used a 2D-Gaussian distribution provided by the 
[#Astropy_Library]_ as the star point spread function as follows:

.. image:: /docs/images/flux_distribution.png   
   :alt: Flux distribution equation
   :align: center

where

.. image:: /docs/images/a_coefficient.png   
   :alt: A coefficient
   :align: center

.. image:: /docs/images/b_coefficient.png   
   :alt: B coefficient
   :align: center

.. image:: /docs/images/c_coefficient.png   
   :alt: C coefficient
   :align: center


:math:`f_p(x,y)` is the star intensity in ADU, C represents the maximum amplitude in ADU, x and y are the coordinates over the image in pixels,
:math:`x_0` and :math:`y_0` are the star coordinates in pixels, :math:`\delta_x` and :math:`\delta_y` are the standard deviation of the Gaussian 
in the x and y directions, respectively; :math:`\theta`; is the rotation angle of the Gaussian.

The standard images created by the simulator have 200 x 200 pixels, the center coordinates of the star was fixed in (:math:`x_0`, :math:`y_0`) = (100,100) pixels; 
the values :math:`\delta_x` = :math:`\delta_y`, and :math:`\theta` = 0. The maximum amplitude C, in ADU, for each image is calculated through

.. image:: /docs/images/photons_flux.png   
   :alt: C coefficient
   :align: center

where :math:`\beta` simulates a constant photon flux over the CCD. Finally, it is created one frame of the image cube by adding thogether the image 
noise and the star flux distribution. So, the image cube is assembled by concatenating several frames created in the previous step. The resulting 
cube is saved in the FITS 16-bit unsigned format.


Running the AIS
===============

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

Prerequisites
-------------

There are some packages that need to be installed before running the software.

* Astropy_
* Photutils_
* Collections_
* JSON_
* Pandas_
* xlrd_
* Scipy_

To install these packages it is suggested to use the pip command as follows :code:`pip install <package_name>`

Installing
----------

Clone this repo using :code:`git clone https://github.com/DBernardes/AIS.git`

Running the tests
-----------------

To run a simple test, you only need to execute the run.py file and the image will be created in your current directory. 
The run.py file will provide to the AIS the basic information for its execution, that is the star flux, in photons/s; the sky flux, 
in photons/pix/s, and the information about the CCD. In particular, the CCD information should be a python dictionary with the control 
parameters used to configure the acquisition of the SPARC4 cameras. These parameters are the Electron Multiplying Mode (em_mode), the 
Electron Multiplying Gain (em_gain), the Pre-amplification (preamp), the Horizontal Shift Speed (hss), the Pixels Binning (bin), and the Exposure 
Time (texp). Below, it is presented the accepted values for each parameter previously described.

- star_flux: greater than zero
- sky_flux: greater than zero
- em_mode: 0 or 1
- em_gain: from 2 to 300
- preamp: 1 or 2
- hss: 0.1, 1, 10, 20, and 30
- bin: 1 or 2
- texp: greater or equal than 1e-5 s


Beyond the paramaters presented before, there are a set of optional paramaters. They are the CCD temperature in celsius degree (ccd_temp), 
the CCD serial number (serial_number), the image bias level in ADU (bias_level), the size of the image cube (cube_size), the number of lines 
and columns of the image in pixels (lines_number and columns_number, respectively), the gain of the CCD in e-/ADU (gain), the dark current noise 
in e- (dark_noise), the read noise in e- (read_noise), the name and directory of the image (image_name and image_dir, respectively), and the 
standard deviation of the gaussian in pixels (gaussian_stddev). If the values of the dark noise, read noise and the gain are not provided, the 
software will set these values based on the operation mode of the CCD.

- ccd_temp: from 0 ºC to -70 ºC
- serial_number: 9914, 9915, 9916, or 9917
- bias_level: integer and greater or equal than 1
- cube_size: integer and greater than 1  
- lines_number: integer and greater or equal than 1
- columns_number: integer and greater or equal than 1
- gain: greater than 0
- dark_noise: greater than zero
- read_noise: greater than zero       
- image_name: string
- image_dir: string
- gaussian_stddev: integer and equal or greater than 1
   

Authors and Contact
====================

* **Denis Bernardes**: 

email: denis.bernardes099@gmail.com 

License
=======

This project is licensed under the MIT License - see the LICENSE_ file for details


References
==========

.. [#SPARC4] Claudia V. Rodrigues, Keith Taylor, Francisco J. Jablonski, Marcelo Assafin, Alex Carciofi, Deonisio Cieslinski, Joaquim E. R. Costa, Ruben Dominguez, Tania P. Dominici, Gabriel A. P. Franco, Damien J. Jones, Antonio Kanaan, René Laporte, Antonio M. Magalhaes, André Milone, José A. Neri, Antonio Pereyra, Luiz A. Reitano, Karleyne M. G. Silva, Cesar Strauss, "Concept of SPARC4: a simultaneous polarimeter and rapid camera in 4 bands," Proc. SPIE 8446, Ground-based and Airborne Instrumentation for Astronomy IV, 844626 (24 September 2012); https://doi.org/10.1117/12.924976

.. [#Bernardes_2018] Bernardes, D. V., Martioli, E., and Rodrigues, C. V., “Characterization of the SPARC4 CCDs”, <i>Publications of the Astronomical Society of the Pacific</i>, vol. 130, no. 991, p. 95002, 2018. doi:10.1088/1538-3873/aacb1e.

.. [#Astropy_Library] The Astropy Collaboration et al 2018 AJ 156 123



.. _Astropy: https://www.astropy.org/
.. _Photutils: https://photutils.readthedocs.io/en/stable/
.. _Collections: https://docs.python.org/3/library/collections.html
.. _JSON: https://www.w3schools.com/python/python_json.asp
.. _Pandas: https://pandas.pydata.org/
.. _xlrd: https://xlrd.readthedocs.io/en/latest/
.. _Scipy: https://www.scipy.org/
.. _LICENSE: https://github.com/DBernardes/AIS/blob/main/LICENSE
.. |Documentation| image:: https://readthedocs.org/projects/ais/badge/?version=latest
	:target: https://ais.readthedocs.io/en/latest/?badge=latest
	:alt: Documentation Status
.. |Unit Tests| image:: https://github.com/DBernardes/AIS/actions/workflows/python-unittests.yml/badge.svg
	:target: https://github.com/DBernardes/AIS/actions/workflows/python-unittests.yml
	:alt: Unit Tests
      
