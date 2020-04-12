Installation
=======================================

Get the pdb_manip_py library from `github`_.

.. code-block:: bash

	git clone https://github.com/samuelmurail/pdb_manip_py.git
	./setup.py install --user

.. _github: https://github.com/samuelmurail/pdb_manip_py_py

Prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. python 3 libraries:
	* `htmd-pdb2pqr`_
	* Numpy
	* Scipy
	* `Os_Command_py`_
	* Sphinx and sphinx-argparse (only for building documentation)
	* Pytest (only for testing purpose)

.. _htmd-pdb2pqr: http://www.poissonboltzmann.org/
.. _Os_Command_py: https://github.com/samuelmurail/os_command_py


Make the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Need `sphinx`_ installed with the argparse sphinx module:

.. code-block:: bash

	pip3 install Sphinx --user
	pip3 install sphinx-argparse --user

You can then build the documentation either in html format or pdf.

.. code-block:: bash

	cd pdb_manip_py/doc
	# For html documentation:
	sphinx-build -b html . _build
	# For pdf documentation:
	sphinx-build -M latexpdf . _build/

.. _sphinx: http://www.sphinx-doc.org

Test installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Launch test with `doctest`_, will check that module’s docstrings are up-to-date by verifying that all interactive examples still work as documented.

.. code-block:: bash

	$ pytest
	================================= test session starts =================================
	platform darwin -- Python 3.6.10, pytest-5.4.1, py-1.8.1, pluggy-0.13.1
	rootdir: /Users/smurail/Documents/Code/pdb_manip_py, inifile: pytest.ini
	collected 33 items

	pdb_manip_py/pdb2pqr.py .                                                       [  3%]
	pdb_manip_py/pdb_manip.py ................................                      [100%]

	================================= 33 passed in 8.25s ==================================
.. _doctest: https://docs.python.org/3/library/doctest.html