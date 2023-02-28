# quantum-matched-filter

For details, see [A quantum algorithm for gravitational wave matched filtering](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.4.023006).

Matched filtering is a technique used in signal processing to detect signals using known template signal. It is used in gravitational wave data analysis to detect gravitational wave signals from astrophysical sources using a template bank of waveforms produced from numerical relativity simulations. Such template banks are often extremely large and comparing each template with the signal can be very computationally expensive if done classicaly. However a quantum approach of shifting through the template bank for matches can offer a significant speed-up. 

`qiskit_example.ipynb` is an example of how quantum matched filtering can be implemented on a quantum computer as a toy model written in the [Qiskit](https://qiskit.org/) software and run on IBM's quantum simulator. `GW150914_example.ipynb` is a step-by-step demonstration of the quantum matched filtering algorithm being used on the first discovered gravitational wave signal GW150914. The data used in this example is pulbically available from [GWOSC](https://www.gw-openscience.org/about/) and the preprocessing steps are carried out in `GW150914_preprocessing.ipynb`.
