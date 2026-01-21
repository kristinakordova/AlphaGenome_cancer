from setuptools import setup, find_packages

setup(
    name="AlphaGenomeAnalysis",  # Or your chosen repository name
    version="0.1.0",
    author="Your Name",
    description="A toolkit for benchmarking AlphaGenome predictions against RNA-Seq data",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "seaborn",
        "scipy",
        "pysam",
        "tqdm",
        "scikit-learn",  # Required for RMSE/MAE calculations
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
