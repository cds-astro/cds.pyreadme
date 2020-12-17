# by G.Landais (CDS) (12/2020)

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gilles.landais", # Replace with your own username
    version="1.3.2",
    author="Gilles Landais",
    description="ReadMe generator package (cdspyreadme)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(exclude=['tests']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_data={'cdspyreadme': ['*.template']},
    include_package_data=True,
    python_requires='>=3.6',
)
