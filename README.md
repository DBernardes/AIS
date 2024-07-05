|Unit Tests| |Documentation| |Code Cov| |Codacy| |Black|

<h1 style="text-align: center">Welcome to the documentation of the Artificial Image Simulator</h1>

<p style="text-align: justify;">
Instituto Nacional de Pesquisas Espaciais (INPE), in collaboration with Laboratório Nacional de Astrofísica (LNA), is developing the Simultaneous Polarimeter and Rapid  Camera in Four Bands (SPARC4), a new astronomical instrument to be installed in the Perkin-Elmer 1.6 m telescope of Picos dos Dias Observatory (in Portuguese, OPD). In this project, we present the Artificial Image Simulator (AIS). AIS was developed to simulate the photometric and polarimetric operating modes of SPARC4, as a function of its spectral response. To obtain this spectral response, we used the theoretical curves provided by the manufacturer of the optical elements and the scientific cameras of the instrument. Also, we considered the spectral response of other systems that are relevant to this characterization, such as the telescope and the atmosphere. The individual optical components that compound the instrument are the elements of the polarimetric module, collimator, dichroic mirrors, and optical focuser cameras. AIS was implemented using the Python programming language. It was used in different activities for the instrument characterization. Among these activities, we can highlight the production of an empirical model to transform the SPARC4 instrumental system into the Sloan Digital Sky Survey (SDSS) system. This model was calibrated using observations made with the instrument. Moreover, the simulator can also be used as an exposure time calculator of SPARC4, for planning scientific observations using the instrument, and investigate the instrumental effects related to polarimetry. 
</p>


![example image](images/example_image.png){width=300 .center}

   

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


.. _LICENSE: https://github.com/DBernardes/AIS/blob/main/LICENSE
.. |Documentation| image:: https://readthedocs.org/projects/ais/badge/?version=latest
	:target: https://ais.readthedocs.io/en/latest/?badge=latest
	:alt: Documentation Status
.. |Unit Tests| image:: https://github.com/DBernardes/AIS/actions/workflows/python-app.yml/badge.svg
	:target: https://github.com/DBernardes/AIS/actions/workflows/python-app.yml
	:alt: Unit Tests
.. |Code Cov| image:: https://codecov.io/gh/DBernardes/AIS/branch/main/graph/badge.svg?token=aPhVaeHkOh
      :target: https://codecov.io/gh/DBernardes/AIS
      
.. |Codacy| image:: https://app.codacy.com/project/badge/Grade/b2724af0f0b043bc84f768659f73fb77    
    :target: https://www.codacy.com/gh/DBernardes/AIS/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DBernardes/AIS&amp;utm_campaign=Badge_Grade

.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

      
