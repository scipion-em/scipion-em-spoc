========================
SPOC plugin
========================

This plugin provide wrappers around several programs of `SPOC <https://github.com/MaximilianBeckers/SPOC>`_.

+------------------+------------------+
| stable: |stable| | devel: | |devel| |
+------------------+------------------+

.. |stable| image:: http://scipion-test.cnb.csic.es:9980/badges/eman2_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/eman2_sdevel.svg


Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version (not avalaible yet)

.. code-block::

    scipion installp -p scipion-em-spoc

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-spoc.git

    * install

    .. code-block::

        scipion installp -p path_to_scipion-em-spoc --devel

SPOC binaries will be installed automatically with the plugin.

* Default installation path assumed is ``software/em/spoc-1.0``, if you want to change it, set *SPOC_HOME* in ``scipion.conf`` file pointing to the folder where the SPOC is installed.

To check the installation, simply run one of the following Scipion tests:

.. code-block::

   scipion3 tests spoc.tests.test_protocols_spoc.TestFscFdrControl

A complete list of tests can also be seen by executing ``scipion test --show --grep emantomo``

Supported versions
------------------

* 1.0

Protocols
---------

* FSC_FDR Control

References
----------

1. Beckers, M., Jakobi, A. J., Sachse, C. (2019) Thresholding of cryo-EM density maps by false discovery rate control. IUCr Journal 6 (1)
