from setuptools import setup, find_packages

setup(
    name="nsforest-cli",
    version="0.1",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        "typer[all]",
        "scanpy",
        "nsforest",
    ],
    entry_points={
        "console_scripts": [
            "nsforest-cli = nsforest_cli.main:app"
        ]
    },
)

