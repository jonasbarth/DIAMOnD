from setuptools import setup, find_packages

setup(
    name='diamond',
    version='0.2',
    packages=find_packages(),
    install_requires=['numpy',
                      'scipy',
                      'pandas',
                      'networkx'],
    entry_points={
            'console_scripts': [
                'diamond = diamond.DIAMOnD:main_function',
            ],
        },
)