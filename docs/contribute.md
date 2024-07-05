# How to Contribute

<div style="background-color: #FFFFE0; padding: 10px; border: 1px solid #E6DB55;">
    <p><strong>Warnings:</strong> this page is a working in progress. Its text is still being modified.</p>
</div>

Thank you for your interest in contributing to the Artificial Image Simulator (AIS) project! This document lists the most common operations you may need to contribute.

## How does the project work?


### Project Structure

The project is divided into three directories: docs, AIS, and tests, where each directory has its specific function.

```
notas_musicais --> escalas.py
notas_musicais --> acordes.py
notas_musicais --> campo_harmonico.py
notas_musicais --> cli.py
```

The CLI code and library are located in notas_musicais. The API documentation for the code is also being developed in notas_musicais using the mkdocstrings tool, following the Google docstring style. Therefore, if you modify anything in the code, remember to update the docstrings as well.

The examples used in the docstrings are also used for testing purposes. Therefore, if you change the output format, remember to update the docstrings.

## About the Library

The entire library uses pure Python without any external library dependencies intentionally to keep the code simple. The function responses are standardized with the return always being a Python dictionary. This standardization facilitates potential future integration into a graphical interface or REST API.

Every time the code is passed between functions during application construction, a pattern has been established to unpack the dictionary in other functions. So, don't worry if you see this code format frequently:

py
Copiar código
notes, degrees = function('arg1', 'arg2').values()

## The CLI
The CLI was built using the Typer library, and you can check its documentation for more details if you want to expand CLI functionalities.

For rich output responses in the application, the Rich library was used. If you want to make changes to the tables generated in the output, you can directly visit the tables documentation page.

The only convention followed regarding the CLI is that a Console object from Rich and a Typer app have already been defined. It would be beneficial if you continued using these objects:

py
Copiar código
from rich.console import Console
from typer import Argument, Typer

...

console = Console()
app = Typer()


## tests
For testing, we are using pytest. Its configurations can be found in the pyproject.toml file at the root of our project.

One important thing to note about the tests is that not all tests are only in the notas_musicais/tests directory. The addopts = "--doctest-modules" flag is being used. Therefore, if you modify something, be aware that docstrings also run tests and serve as the basis for API documentation. So, be careful with changes.

Test coverage is automatically generated with pytest-cov and is displayed when the test task is executed:

bash
Copiar código
task tests
Similarly, linters are prerequisites for these tests.

## Documentation
The entire documentation is based on the use of mkdocs with the mkdocs-material theme.

mermaid
Copiar código
flowchart
    . --> docs
    . --> mkdocs.yml
	docs --> arquivos.md
	docs --> api
	docs --> assets
	docs --> templates
	docs --> stylesheets
All configuration can be found in the mkdocs.yml file at the root of the repository.

Various techniques are also being used to complement the documentation, such as templates from Jinja where instructions may repeat. If you encounter blocks like:

html
Copiar código
{ %  % }
You'll know it's a template.

The templates are defined in the /docs/templates directory. In some cases, however, they may be called by variables with command.run that appear in almost all documentation files. These macros are created with mkdocs-macros and are defined in the mkdocs configuration file:

yaml
Copiar código
extra:
  commands:
    run: poetry run notas-musicais

## API Documentation
The API documentation is being done within the code modules. Therefore, files in the docs/api directory have a tag:

md
Copiar código
This means that the code contained in the docstrings will be used in this block. The mkdocstrings plugin is used to handle this.

The module documentation follows the Google docstring format, which is the library's standard.

## Tools
This project basically uses two tools as its foundation:

Poetry: For environment management and library installation
Taskipy: For automation of routine tasks, such as running tests, linters, documentation, etc.
So, ensure that Poetry is installed for this contribution:

bash
Copiar código
pipx install poetry

## Steps to Execute Specific Tasks
Here are commands you can use to perform routine tasks, such as cloning the repository, installing dependencies, running tests, etc.

## How to Clone the Repository
bash
Copiar código
git clone https://github.com/dunossauro/notas-musicais.git

## How to Install Dependencies
bash
Copiar código
poetry install

## How to Run the CLI
bash
Copiar código
poetry run notas-musicais [subcommand]

## How to Run Code Verification
bash
Copiar código
task lint
How to Run Tests
bash
Copiar código
task test

## How to Run Documentation
bash
Copiar código
task docs
Tasks You Can Contribute To
{% include "templates/todos.md" %}

## Didn't find what you need?
If you didn't find what you need, you can open an issue in the project to report what you can't do or what needs better documentation.

## Continuous Improvement
This document can be improved by anyone interested in enhancing it. So, feel free to provide more tips for people who want to contribute too! :heart: