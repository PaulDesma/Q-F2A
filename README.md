Q+F2A
==========

Python code used to compute the frequency dependent admittance of Electric Double Layer using Brownian Dynamics using forces on ions and, inspired from Q2Z of Giovanni Pireddu (https://github.com/gpireddu/Q2Z).

It takes the time series of the charge and sum of the forces over each ion species to compute the complexe admittance spectra. It is based on the linear response theory and relies on Laplace transform using Filon integrations.
 
 ---
# Reference

[Paul Desmarchelier and Benjamin Rotenberg, arXiv, 2026](https://arxiv.org/abs/2607.00932)

Bibtex:
```
@misc{desmarchelier2026dynamicschargefluctuationsnanocapacitors,
      title={Dynamics of charge fluctuations in nanocapacitors: effects of salt concentration and electrode metallicity from Brownian dynamics}, 
      author={Paul Desmarchelier and Benjamin Rotenberg},
      year={2026},
      eprint={2607.00932},
      archivePrefix={arXiv},
      primaryClass={physics.chem-ph},
      url={https://arxiv.org/abs/2607.00932}, 
}
```


---
# Dependencies
Q+F2A requires the following packages to run:

* numpy
* scipy
* pyfilon.pyfilon (https://pypi.org/project/pyfilon/, https://github.com/alexhroom/pyfilon)

# Usage
The electrode charge and sum of the forces time series should be in the same folder and named ```total_charges.out.gz``` and ``` ionic_current.out.gz```. An example is provided in ```Example/``` for testing purposes only.

Simply run:

```python3 MainAdmittance.py```

The complex admittance and associated variance will be outputed in a file named ``` Admlen$[correlationLength]nfreq$[numberOfFrequencysampled]```  and ``` VarAdmlen$[correlationLength]nfreq$[numberOfFrequencysampled].data``` 
# Notes

* The code assumes a MetalWalls output format [https://gitlab.com/ampere2/metalwalls/-/wikis/output-files#total_charges.out](https://gitlab.com/ampere2/metalwalls/-/wikis/output-files#total_charges.out)
* Parameters such as the temperature, duration of time steps and the list of desired frequencies are hard-coded, for the moment.
* The code is not tested for Python versions other than Python3


# Files
- funcAdmittance.py
 	  - function called in mains linked to data post process
- MainAdmittance.py
  	 - post process of the charge /ion fluctuation 
