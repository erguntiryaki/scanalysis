from setuptools import setup

setup(
    name='schelper',
    version='0.1.0',
    packages=[''],
    package_dir={'': 'helper'},
    url='',
    license='',
    author='ergun',
    author_email='erguntiryaki27@gmail.com',
    description='Helper functions for scRNA-seq analyses',
    install_requires=['scanpy', 'pandas', 'seaborn']
)
