[![Documentation Status](https://readthedocs.org/projects/ais/badge/?version=latest)](https://ais.readthedocs.io/en/latest/?badge=latest)
[![Unit tests](https://github.com/DBernardes/AIS/actions/workflows/python-app.yml/badge.svg)](https://github.com/DBernardes/AIS/actions/workflows/python-app.yml)
[![Code cov](https://codecov.io/gh/DBernardes/AIS/branch/main/graph/badge.svg?token=aPhVaeHkOh)](https://codecov.io/gh/DBernardes/AIS)
[![Codacy](https://app.codacy.com/project/badge/Grade/b2724af0f0b043bc84f768659f73fb77)](https://www.codacy.com/gh/DBernardes/AIS/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DBernardes/AIS&amp;utm_campaign=Badge_Grade)
[![DOI](https://zenodo.org/badge/346361531.svg)](https://zenodo.org/doi/10.5281/zenodo.12808262)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPi license](https://badgen.net/pypi/license/pip/)](https://pypi.org/project/pip/)
 [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)

<h1 style="text-align: center">Artificial Image Simulator of SPARC4</h1>


The Simultaneous Polarimeter and Rapid Camera in Four Bands ([SPARC4](https://coast.lna.br/home/sparc4)) is a new astronomical instrument developed by *Instituto Nacional de Pesquisas Espaciais* (INPE), in collaboration with *Laboratório Nacional de Astrofísica* (LNA), that is currently installed on the Perkin-Elmer 1.6 m telescope of Picos dos Dias Observatory (in Portuguese, OPD). SPARC4 was developed to allow photometric and polarimetric acquisitions in the g, r, i and, z bands of the Sloan Digital Sky Survey ([SDSS](https://www.sdss.org/)). The acquisition is done using four Framre Transfer (FT) Electron Multplying CCD (EMCCD) cameras produced by Oxford Instruments, one for each optical band of the instrument. 

In this project, we present the Artificial Image Simulator ([AIS](https://zenodo.org/doi/10.5281/zenodo.12808262)), a software developed using the Python programming language for the creation of artificial images, similar to those images acquired using SPARC4. The process for the creation of these images can be split into the modeling of the source Spectral Energy Distribution (SED), the application of the spectral response of optical systems considered in the light path, and the creation of the images themselves. 

For modeling the source SED, the user can choose between using a black-body distribution or a library of spectrophotometric standard stars. Besides, the user can provide their own SED if necesssary. Then, AIS provides the option of applying a polarization to the source SED. This polarization can be a linear, circular, or elliptcal flat distributions, or a distribution given by the Serkowski's polarimetric curve. 

Once the source SED is obtained, the spectral response of the optical systems considered in the source's light path are applied. These systems are the atmosphere, the telescope and the instrument itself. This procedure results in the spectral distribution of photoelectrons that reach the EMCCD cameras. Then, the number of photoelectrons per second acquired by these cameras can be found by integrating this distribution as a function of the wavelength and also considering the area of the telescope primary mirror. 

The artificial images are composed of three another images: the background image, the source image, and the source noise distribution image. The background image is composed of the photons originated from the sky, the background level given by the camera, as well as the noise distribution related to these components. The source image is composed of the Point Spread Function (PSF) of the source. This PSF is modeled as a Gaussian 2D distribution, using the number of photoelectrons per second obtained in the previous step. The source noise distribution image is obtained by multiplying pixel-by-pixel an image given by the Poisson distribution of the source photons and the Gaussian distribution obtained for the source PSF. Finally, the articial image created by AIS is obtained by summing the images created for the background, the source and the source noise distribution.

## Installation

Use the `git clone` command to install the repository into your local computer.

```bash
git clone https://github.com/DBernardes/AIS.git
```

Then, use the following commands to install AIS into your local Python environment.

```bash
cd AIS
pip install -e .
```

## Usage

Once you have this project installed, you will be able to run the file `example.py` found into the `./scritps` folder. For that, use:

```bash
cd scripts
python -m example.py
```
This script will create a basic artificial image into the `./scripts/FITS` folder.


## Contributing

If you want to contribute to this project, please take a look to the following guidelines...

## Authors and acknowledgment

- **Denis Bernardes (main developer)**
  - [![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=flat&logo=github&logoColor=white)](https://github.com/DBernardes) 
[![LinkedIn](https://img.shields.io/badge/linkedin-%230077B5.svg?style=flat&logo=linkedin&logoColor=white)](www.linkedin.com/in/denisbernardes) 
[![Gmail](https://img.shields.io/badge/Gmail-D14836?style=flat&logo=gmail&logoColor=white)](mailto:denis.bernardes099@gmail.com)
  - Affiliation: Instituto Nacional de Pesquisas Espaciais (INPE).
  - Address: 1758 Astronautas Avenue, Jardim da Granja, São José dos Campos, São Paulo, Brazil.

- **Claudia Vilega Rodrigues (supervisor)**
  - Affiliation: Instituto Nacional de Pesquisas Espaciais (INPE).
  - Address: 1758 Astronautas Avenue, Jardim da Granja, São José dos Campos, São Paulo, Brazil.


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details




      
