Fabber models for quantitative BOLD MRI
=======================================

These models use the Fabber_
Bayesian model fitting framework [1]_ to implement a model
for quantitative BOLD MRI.

.. note::
    The Quantiphyse_ application contains a widget for performing
    qBOLD analysis which uses ``FABBER_QBOLD`` internally.
    
Getting FABBER_QBOLD
--------------------

The qBOLD models are included in FSL_. We
stongly recommend version 6.0.1 or later.

If you need an updated version of the model which has not yet been released to
FSL, you will either need to 
`build from source <https://fabber-core.readthedocs.io/en/latest/building.html#building-new-or-updated-model-libraries>`_ 
using an existing FSL 6.0.1 or later installation, or download 
the pre-built `Fabber bundle <https://fabber-core.readthedocs.io/en/latest/getting.html#standalone-fabber-distribution>`_ 
which contains the latest qBOLD release alongside other models in a standalone package.

The qBOLD model
---------------

Examples
--------

References
----------

.. [1] *Chappell, M.A., Groves, A.R., Woolrich, M.W., "Variational Bayesian
   inference for a non-linear forward model", IEEE Trans. Sig. Proc., 2009,
   57(1), 223–236.*

.. _Fabber: https://fabber-core.readthedocs.io/

.. _FSL: https://fsl.fmrib.ox.ac.uk/fsl/

.. _Quantiphyse: https://quantiphyse.readthedocs.io/
