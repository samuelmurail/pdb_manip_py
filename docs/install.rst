Installation
=======================================

Using Pypi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip3 install pdb_manip_py

Using Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   conda install pdb_manip_py -c conda-forge


From source code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The sources for Docking Python can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/samuelmurail/pdb_manip_py

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/samuelmurail/pdb_manip_py/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install

.. _Github repo: https://github.com/samuelmurail/pdb_manip_py
.. _tarball: https://github.com/samuelmurail/pdb_manip_py/tarball/master


Test installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Launch test with `doctest`_, will check that module’s docstrings are up-to-date by verifying that all interactive examples still work as documented.

.. code-block:: bash

	$ pytest
	=============================== test session starts ================================
	platform darwin -- Python 3.8.2, pytest-5.4.1, py-1.8.1, pluggy-0.13.1
	rootdir: /Users/smurail/Documents/Code/pdb_manip_py, inifile: pytest.ini
	plugins: cov-2.8.1
	collected 40 items

	pdb_manip_py/pdb2pqr.py .                                                    [  2%]
	pdb_manip_py/pdb_manip.py .......................................            [100%]

	================================= warnings summary =================================
	pdb_manip_py/pdb_manip.py::pdb_manip_py.pdb_manip.Coor.align_principal_axis
	  /Users/smurail/opt/miniconda3/envs/docking/lib/python3.8/site-packages/scipy/spatial/transform/rotation.py:1953: UserWarning: Optimal rotation is not uniquely or poorly defined for the given sets of vectors.
	    warnings.warn("Optimal rotation is not uniquely or poorly defined "

	-- Docs: https://docs.pytest.org/en/latest/warnings.html
	========================== 40 passed, 1 warning in 11.64s ==========================

.. _doctest: https://docs.python.org/3/library/doctest.html