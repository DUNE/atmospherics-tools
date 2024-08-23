# Atmospherics-tools

This repository aims at regrouping useful analysis tools in the context of the various studies made on DUNE atmospheric neutrinos. The tools are spread among several folders that each have their README.

## [anatree](anatree)

Contains a python class to read-in the anatrees with Polars. Allows to separate the various information in different sub-dataframes: nu, geant, reco, pfps.
Also contains some generic python analysis tools.

## [example-notebooks](example-notebooks)
Contains several example notebooks that demonstrate some of the possible analysis using files from the atmospherics production and some tools of this repo.

## [fast-osc-feedback](fast-osc-feedback)
Attempt at developing a simple python tool to be able to get fast feedback on the impact of the reconstruction quality on the DUNE sensitivity to the oscillation parameters. Still WIP

## [weight_calculation](weight_calculation)
Tools to perform the event weight calculation on the atmospherics samples. Allows to compute different flux weights based on the incoming nu type as well as oscillated weights taking into account the transition probabilities from one flavour to the other

## [xsec_systs_calculation](xsec_systs_calculation)
Tool based on nusystematics allowing to compute the impact of xsec parameters shifts on event weights. This is used to compute the xsec systs splines that are inputed into MaCh3 for the oscillation analysis

## General comments

Please, feel free to make PR in order to add more tools, more examples, more documentation to this repository!
