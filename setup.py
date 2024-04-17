import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="crispat",
    version="0.1.0",
    author="Jana M. Braunger",
    author_email="jana.braunger@cos.uni-heidelberg.de",
    description="crispat: CRISPR guide Assignment Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gtca/mofax",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=[],
)