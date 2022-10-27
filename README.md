# METAS UncLib Python Wrapper

For METAS UncLib see

- <https://www.metas.ch/unclib>
- <https://pypi.org/project/metas-unclib/>

## Tutorial Notebook

- [METAS UncLib Python Tutorial](metas_unclib/jupyter_notebooks/metas_unclib_python_tutorial.ipynb)

## Example Notebooks

- [Right Triangle Example](metas_unclib/jupyter_notebooks/right_triangle_example.ipynb)
- [Resistor Cube Example](metas_unclib/jupyter_notebooks/resistor_cube_example.ipynb)

## Wrapper

- [METAS UncLib Pyhton Wrapper](metas_unclib/metas_unclib.py)

## Developer workflow
- Install [METAS UncLib](https://www.metas.ch/metas/en/home/fabe/hochfrequenz/unclib.html)
    - Delete PYTHONPATH environment variable to enable development within virtual environment
- Fork repo on Github
- Clone fork onto local machine and navigate to folder
- Create a virtual environment, activate and update pip
- Copy UncLib DLLs into package folder `metas_unclib`
    - Program Files (x86)/METAS/Metas.UncLib/\*.dll*
    - Program Files (x86)/METAS/Metas.UncLib/mkl_custom/**
- pip install -e .[dev]
- Configure editor to use `black`, `pylint` and `pytest`
- Create a feature branch
- Add feature and tests
- Run pytest, black and pylint
- Build package: python -m build
- Test wheel in separate virtual environment
- Push feature branch to fork and create a PR
