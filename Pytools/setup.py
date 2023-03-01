import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MattTools",
    version="0.1",
    author="Matthew Muller",
    author_email="matt.alex.muller@gmail.com",
    description="Some persional functions",
    url="https://github.com/mattmuller0/scripts/tree/main/Pytools",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=[
        "matplotlib >= 3.3.4",
        "numpy >= 1.20.1",
        "pandas >= 1.2.2",
        "scipy >= 1.6.0",
        "seaborn >= 0.11.1",
        "sklearn >= 1.1.3",
    ]
)