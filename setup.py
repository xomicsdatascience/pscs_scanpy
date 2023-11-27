import setuptools
import os


def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()


def get_version():
    """Gets the version from __init__ without importing it"""
    here_dir = os.path.dirname(__file__)
    pscs_init = os.path.join(here_dir, 'pscs_scanpy', '__init__.py')
    init_contents = readfile(pscs_init).split('\n')
    version = "unknwon"
    for line in init_contents:
        if "__version__" in line:
            version = line.split('=')[1].replace('"', '').strip()
            break
    return version


def get_requirements():
    """Loads the contents of requirements.txt and returns them as a list."""
    here_dir = os.path.dirname(__file__)
    req_path = os.path.join(here_dir, 'requirements.txt')
    if os.path.exists(req_path):
        return readfile(req_path).splitlines()
    else:
        return None


setuptools.setup(
    name="pscs_scanpy",
    version=get_version(),
    author="Alexandre Hutton",
    author_email="",
    description="Scanpy port for PSCS",
    long_description="",
    url="",
    classifiers=[
        "Development Status :: 2 - Alpha",
        'Intended Audience :: Science/Research',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: AGPL 3"
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    python_requires='>=3.11',
    license="AGPL 3",
    install_requires=get_requirements(),
)
