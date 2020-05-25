import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="msp-IrTrez",
    version="0.0.1",
    author="Tristan Dijkstra, Andrea Battegazzore, Hugo Chassagnette",
    author_email="tristan.dijkstra@gmail.com",
    description="MSP: two-body astronautical model written in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/IrTrez/MSP",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL v3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
