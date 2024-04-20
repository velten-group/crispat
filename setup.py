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
    url="https://github.com/velten-group/crispat/tree/main",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=["anndata>=0.8.0", "matplotlib>=3.7.1", "numpy>=1.24.3", "pandas>=1.5.3", "pyro_ppl>=1.8.6", 
                      "scanpy>=1.9.3", "scipy==1.9.3", "seaborn==0.12.2", "setuptools==69.5.1", "torch==2.0.1",
                      "tqdm==4.65.0"],
)