"""Setup file of the Artificial Images Simulator of the SPARC4 instrument."""

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name="AIS",  # Replace with your own username
    version="1.0",
    license="MIT",
    author="Denis Varise Bernardes",
    author_email="denis.bernardes099@gmail.com",
    description="""This is the Artificial Images Simualtor (AIS) developed to
    create artificial images of stars, simulating those images acquired by the
    SPARC4 CCD cameras in astronomical observations""",
    long_description=long_description,
    long_description_content_type="restructured text",
    url="https://github.com/DBernardes/AIS.git",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.10",
)
