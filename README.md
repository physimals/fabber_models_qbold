# qBOLD Forward Model

An implementation of the qBOLD model (based on He & Yablonskiy, 2007, and others) for FABBER.

Everything necessary should be in `fwdmodel_qbold.cc`, although there is probably a more elegant, object-oriented method to dealing with the different model options (currently there are just a bunch of boolean inputs that the user can specify when running it). 

The `buildfabber` bash script should compile everthing (if `fabber_core` is located correctly and `CMakeLists.txt` is right), and will produce an executable called `fqbold` which should then be called with the `runfabber` bash script. This is designed to take a `.fab` file called `runX.fab` with all the input arguments in one place (so that you don't have to type them in the command line)

# the qBOLD Model

The qBOLD model (as implemented here) is based on asymmetric spin-echo data, which has parameters TE (echo time) and tau (spin echo offset). Each volume of the input NIFTI file needs a tau and a TE (or one TE for the whole set) which should be listed in `runX.fab` 

The model describes ASE signal originating from up to three biological compartments: tissue (grey matter), blood (intravascular) and CSF (or extracellular fluid). By default, only the tissue compartment is evaluted. The user can add boolean flags `include_intra` (or `motional_narrowing`) to evaluate the blood compartment using either the powder model (Sukstanskii & Yablonskiy, 2001 - default) or the motional narrowing model (Berman & Pike, 2017). The user can also specify `include_csf` to evaluate the CSF compartment.

Another option is `ignore_T1`, if this is on, then the total signal will be the sum of the signal from each compartment, weighted by that compartment's volume fraction (DBV is the blood compartment volume fraction, and lambda is the CSF compartment volume fraction). If T1 is not being ignored, the compartments are weighted by their steady state magnetization. In this case, the user must also supply TR (repetition time) and TI (inversion time) parameters to calculate magnetization. These must be the same for all volumes in the input file.

The key parameters of the model are OEF and DBV, although inferring on both of these simultaneously is unreliable. It is better to infer on R2p and DBV and then calculate OEF from these later. S0 (the baseline signal) must always be inferred on. If you have data that contains a CSF compartment, you can also infer on lam (lambda) and/or dF (CSF frequency shift), but these won't work well unless you supply some other information. If your data has multiple TEs, you can infer on R2t (which is 1/T2 of tissue) and/or R2e (1/\T2 of CSF), although I haven't properly tested this. Also, you can infer on Hct (fractional hematocrit) but I don't recommend doing so since it will also be a confound with DBV and R2p. 