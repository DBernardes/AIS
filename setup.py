"""Setup file of the Artificial Images Simulator of the SPARC4 instrument."""

import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name="Artificial-Images-Simulator",  # Replace with your own username
    version="1.0",
    license="MIT",
    author="Denis Varise Bernardes",
    author_email="denis.bernardes099@gmail.com",
    description="""This is the Artificial Images Generator (AIG) developed to 
    create artificial star images, simulating the star images acquired by the 
    SPARC4 CCD cameras in astronomical observations""",
    long_description=long_description,
    long_description_content_type="markdown",
    url="https://github.com/DBernardes/AIS.git",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.7",
)
