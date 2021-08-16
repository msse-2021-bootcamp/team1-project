import setuptools

setuptools.setup(
    name="mcsim",
    version="0.0.1",
    author="Anna Weber, Joyce Kim",
    author_email="annaweberdev@gmail.com",
    description="A package for running a Monte Carlo simulation",
    long_description="A Python package which performs MC simulation",
    long_description_content_type="text/markdown",
    url="https://github.com/USERNAME/REPONAME",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.9",
)