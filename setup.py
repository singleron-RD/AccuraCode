import setuptools
from accuracode.__init__ import __VERSION__, ASSAY_DICT

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fp:
    install_requires = fp.read()

entrys = ['accuracode=accuracode.accuracode:main',]
for assay in ASSAY_DICT:
    entrys.append(f'multi_{assay}=accuracode.{assay}.multi_{assay}:main')
entry_dict = {
        'console_scripts': entrys,
}


setuptools.setup(
    name="accuracode",
    version=__VERSION__,
    author="liji",
    author_email="liji@singleronbio.com",
    description="AccuraCode Analysis Pipelines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/singleron-RD/AccuraCode",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=install_requires,
)
