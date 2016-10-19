from setuptools import setup

setup(
    name='pvpy',
    version='1.0.alpha',
    classifiers=[
        'Development Status :: 1 - Development',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific Computing :: Solar',
      ],
    description='Some routines relating to the ASTM G173-03 solar spectrum',
    url='https://github.com/soamaven/pvpy.git',
    author='Colton Bukowsky',
    author_email='cb@caltech.edu',
    license='MIT',
    packages=[
        'pvpy',
    ],
    install_requires=[
        'numpy',
        'scipy'
    ],
    include_package_data=True,
    zip_safe=False
)
