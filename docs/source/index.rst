QOCCM
========

QOCCM (Quick Ocean Carbon Cycle Model) is a powerful one-dimensional ocean carbon cycle model. QOCCM is accurate -- models just like QOCCM are used to emulate the behavior of complex Earth System Models. And QOCCM is quick -- generating hundreds of years of ocean carbon cycle simulation in seconds on a laptop.

The QOCCM package includes code for running idealized experiments such as ocean uptake with fixed buffer capacity. See :doc:`Methodology <methodology>` for details. Please cite Ridge and McKinley (2020) if you use these idealized experiments in your own research.


Installation
------------

Install QOCCM by running:

.. code-block:: bash

    pip install qoccm

If you want to be able to run the :doc:`examples`, install from GitHub

.. code-block:: bash

	git clone https://github.com/sridge/qoccm.git
	cd qoccm
	python setup.py install qoccm

Contents
--------

.. toctree::
   methodology
   examples
   
.. automodule:: qoccm
   :members:

Contribute
----------

- Issue Tracker: `QOCCM GitHub issue tracker <https://github.com/sridge/qoccm/issues>`_ 
- Source Code: `QOCCM GitHub repository <https://github.com/sridge/qoccm>`_ 

Support
-------

If you are having issues, please let me know: sridge@ldeo.columbia.edu

License
-------

The project is licensed under the MIT license.
