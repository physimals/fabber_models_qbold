# qBOLD Forward Model

An implementation of the qBOLD model (based on He & Yablonskiy, 2007, and others) for FABBER.

Everything necessary should be in `fwdmodel_qbold.cc`, although there is probably a more elegant, object-oriented method to dealing with the different model options (currently there are just a bunch of boolean inputs that the user can specify when running it). 

To build with an existing FSL installation see the generic fabber build instructions at:

    https://fabber-core.readthedocs.io/en/latest/building.html

# the qBOLD Model

The qBOLD model (as implemented here) is based on asymmetric spin-echo data, which has parameters TE (echo time) and tau (spin echo offset). Each volume of the input NIFTI file needs a tau and a TE (or one TE for the whole set) which is set for example by ``--tau1=0.5 --tau2=0.6...`` A single TE can be set for example using ``--te=0.0082`` or one for each volume, e.g. ``--te1=0.0082 --te2=0.009...``.

The model describes ASE signal originating from up to three biological compartments: tissue (grey matter), blood (intravascular) and CSF (or extracellular fluid). By default, only the tissue compartment is evaluted. The user can add boolean flags `--incintra` (or `--motion-narrowing`) to evaluate the blood compartment using either the powder model (Sukstanskii & Yablonskiy, 2001 - default) or the motional narrowing model (Berman & Pike, 2017). The user can also specify `--inccsf` to evaluate the CSF compartment.

Another option is `ignore-t1`, if this is on, then the total signal will be the sum of the signal from each compartment, weighted by that compartment's volume fraction (DBV is the blood compartment volume fraction, and lambda is the CSF compartment volume fraction). If T1 is not being ignored, the compartments are weighted by their steady state magnetization. In this case, the user must also supply the repetition time TR (e.g. ``--tr=3.0``) and the inversion time TI (e.g. ``--ti=1.2``) to calculate magnetization. These must be the same for all volumes in the input file.

The key parameters of the model are the oxygen extraction fraction ``oef`` and the deoxygenated blood volume ``dbv``, although inferring on both of these simultaneously is unreliable. It is better to infer on the modified T2 rate ``r2p`` and ``dbv`` and then calculate OEF from these later. ``sig0`` (the baseline signal) should always be inferred on. 

If you have data that contains a CSF compartment, you can also infer on ``lam`` (lambda) and/or ``df`` (CSF frequency shift), but these won't work well unless you supply some other information. If your data has multiple TEs, you can infer on ``r2t`` (which is 1/T2 of tissue) and/or ``r2e`` (1/T2 of CSF), although I haven't properly tested this. Also, you can infer on ``hct`` (fractional hematocrit) but I don't recommend doing so since it will also be a confound with ``dbv`` and ``r2p``.
