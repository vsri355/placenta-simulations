==========
The theory
==========

The theory behind this model is based on the integration of models described in detail in:

 - `Clark et al. Interface Focus, 2015 <http://rsfs.royalsocietypublishing.org/content/5/2/20140078>`_.
 - `Erlich et al. Interface Focus 2019 <https://royalsocietypublishing.org/doi/full/10.1098/rsfs.2019.0021>`_.
 - `Erlich et al. Sci Adv 2019 <https://www.science.org/doi/10.1126/sciadv.aav6326>`_.
 - `Byrne et al. J Theoret  Biol, 2021 <https://www.sciencedirect.com/science/article/pii/S0022519321000527>`_.

Please remember to cite these papers if using our softwares.

The basics
----------

The model aims to connect macro-scale perfusion simulations (e.g. Clark et al 2015) with micro-scale models that define
geometric determinants of terminal villous exchange capacity (e.g. Erlich et al 2019a/b).

The theory behind the perfusion model is detailed in :ref:`steady`. However, in this model there are some
modifications:

- First, 'capillary convolutes' as defined by Clark et al. are replaced by :math:`n_S` series and :math:`n_P` parallel terminal villi units as defined by Erlich et al. These terminal villous units have a resistance associated with them which reflects their geometric structure.
- Second, vessels were assumed to be compliant: :math:`r=r_0(1+3r_0 P_{tm}/4Eh_0)` (where r is radius, h is wall thickness, E is elastic modulus and :math:`P_{tm}` is transmural pressure.
- Finally a non-linear rheology model is employed (Pries, A., TW, S., Gaehtgens, P., and JF, G. (1990). Blood flow in microvascular networks: experiments and simulations. Circ Res, 67:826–834)

Nutrient exchange is assumed to occur only in the terminal villi, and exchange modelling follows  `Erlich et al. Sci Adv 2019 <https://www.science.org/doi/10.1126/sciadv.aav6326>`_ directly. This model provides
an exchange flux and nutrient concentration supplied to each terminal villous which is then
summed over the anatomical tree network to provide a total nutrient delivery/update to the fetal circulation via the placenta.




