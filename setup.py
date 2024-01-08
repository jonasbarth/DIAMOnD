from setuptools import setup, find_packages

setup(
    name='diamond',
    version='0.3',
    packages=find_packages(),
    install_requires=['numpy',
                      'scipy',
                      'pandas',
                      'tqdm',
                      'ndex2'
                      'networkx'],
    entry_points={
            'console_scripts': [
                'diamond = diamond.DIAMOnD:main_function',
            ],
        },
)