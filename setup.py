from setuptools import setup, find_packages

setup(
    name='diamond',
    version='0.1',
    packages=find_packages(),
    install_requires=['numpy',
                      'scipy',
                      'networkx'],
    entry_points={
            'console_scripts': [
                'diamond = diamond.DIAMOnD:main_function',
            ],
        },
)