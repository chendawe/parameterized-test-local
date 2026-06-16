- __Environment Dependencies__
    - Install the package dependencies:
    ```sh
    # conda env export > environment.yml
    conda env create -f full_env.yml
    # pip freeze > requirements.txt
    pip install -r requirements.txt
    ```
- __Codes Roadmap__
    - `WF4Py` codes:
        - The code for implementing the PyCBC plugin: run the code in `install_plugin.txt` in [Parameterized_PN_Tests/pycbc-WF4Py_dchi_modified_IMRPhenomXAS](Parameterized_PN_Tests/pycbc-WF4Py_dchi_modified_IMRPhenomXAS)
        - `PyCBC` modifications based on the `IMRPhenomXAS` implementation of `WF4Py`:
            - The `IMRPhenomXAS.py` code: [pycbc/IMRPhenomXAS.py](pycbc/IMRPhenomXAS.py)
            - The `waveform.py` code (for generating the corresponding time domain waveform): [pycbc/waveform.py](pycbc/waveform.py)
            - USAGES: Directly replace the corresponding file in WF4Py and PyCBC.
    - (_DEPRECATED_) [LALSuite codes](./preload_c_library): codes for implementing the non-gr PN coefficients for `LALSuite` by `preload`.
    - For the codes for running on the server for parameter estimation, please see [server_backup.tar.gz](server_backup.tar.gz).


- __Data Analysis Pipeline__
    - Notebooks for reading `.hdf`s, visualizing and outputting `.csv`s for our modified `IMRPhenomXAS`'s PE results on O1/O2 ten GW events: [Parameterized_PN_Tests/show_PE_results](Parameterized_PN_Tests/show_PE_results)
        - csv path: [Parameterized_PN_Tests/_csv_my_deviations](Parameterized_PN_Tests/_csv_my_deviations)
        - `Mathematica` notebooks for calculating the critical probabilities for each event: [Parameterized_PN_Tests/mathematica](Parameterized_PN_Tests/mathematica)
    - codes for the same functionality based on the LVK official data release for the O1/O2 five GW events: [Parameterized_PN_Tests/gwtc-1_official](Parameterized_PN_Tests/gwtc-1_official)
        - csv path: [Parameterized_PN_Tests/csv](Parameterized_PN_Tests/csv)
        - `Mathematica` notebooks for calculating the critical probabilities for each event: [Parameterized_PN_Tests/mathematica/GWTC-1/multi](Parameterized_PN_Tests/mathematica/GWTC-1/multi)

    - Utils:
        - `Mathematica` notebook for PN coefficient computations, inverse-mapping implementations and plottings: [Public/mathematica/PN_inverse/InversePNUtils.nb](Public/mathematica/PN_inverse/InversePNUtils.nb)