# RioBermejoPOC
###### This repository contains R scripts and field data used to model fluvial POC transit and turnover for the Rio Bermejo, Argentina.
This project implements SoilR to model the transfer of particulate organic carbon (POC) from upstream to downstream through a large alluvial river system.
We model the river POC flux and fit the model output delta14C and OC load to field measurements of these values.
###### Files:
'ber_Cdat_indiv.csv' contains 'carbon stocks' for individual suspended sediment samples collected along the Rio Bermejo.
'ber_Cdat_DI.csv' contains depth-integrated 'carbon stocks' for samples collected along the Rio Bermejo.
'ber_C14dat_indiv.csv' contains measured delta14C values for individual samples collected along the Rio Bermejo.
'ber_C14dat_DI.csv' contains depth-integrated delta14C values for samples collected along the Rio Bermejo.
'marisa_bermejo_5pool_17pars_DI.R' is an R script containing code that parameterizes a carbon cycle model to fit the field data.

'SeriesPools14.R' is an R script containing code that models carbon and radiocarbon over time.

###### Variables:
- *dist* = distance downstream along the channel (km)
- *time* = sediment transit time (yr)
- *C_kg* = suspended POC load (kg)
- *C_kg_SD* = uncertainty on the suspended POC load (kg)
- *d14C* =  measured delta14C values
- *d14C_SD* = standard deviation of measured delta14C values (per mille)
