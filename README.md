# quantum-matched-filter
A quantum matched filter for all of your gravitational wave needs!


Matched filtering is a technique used in signal processing to correlate a known template signal to an observed unkown signal. It is used in gravitational wave data analysis to detect gravitational wave signals from astrophysical sources using a template bank of waveforms produced from numerical relativity simulations. Such template banks are often extremely large and correlating each template with the signal can be very computationally expensive if done classicaly. However a quantum approach of shifting through the template bank for matches could offer a significant speed-up. 

`QMF.ipynb` details how Grover's algorithm works for finding templates, and demonstrates this with an example of performing template matching with sine wave templates. Similarly, `single_qubit.ipynb` details the single qubit approach (currently doesn't work but we may fix later). The repository also contains a report detailing the computational cost of the quantum algorithm in comparison to the classical case.
