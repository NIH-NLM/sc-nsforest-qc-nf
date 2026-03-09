from setuptools import setup, find_packages

setup(
    name="nsforest-cli",
    version="1.0.0",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        "typer[all]",
        "scanpy",
        "numpy",
        "pandas",
        "scipy",
        "plotly",
        "kaleido",  # For plotly image export
    ],
    entry_points={
        "console_scripts": [
            "nsforest-cli=nsforest_cli.main:app"
        ]
    },
)
