# MF2005-SDA-HDF

## A MODFLOW package to linearize Stream Depletion Analysis ##
The conventional numerical method is computationally intensive and prone to numerical noises for stream depletion analyses using MODFLOW. In this study, a new MODFLOW package has been developed to improve the computational efficiency and reduce the noises for each simulation. Using the assumption of unchanged flow coefficients between the baseline and scenario runs, the nonlinear groundwater flow system is linearized for solving the flow equations. The new package has been successfully applied to a regional groundwater model in Nebraska. Our simulation experiments show that the numerical noises, commonly identified in conventional approach, have been significantly reduced and a significant speedup has been achieved for regional groundwater models.

Please cite these two papers if you use SDA in your work.

Ou, G., Li, R., Pun, M., Osborn, C., Bradley, J., Schneider, J., Chen, X.-H., 2016. A MODFLOW package to linearize stream depletion analysis. Journal of Hydrology 532, 9–15. https://doi.org/10.1016/j.jhydrol.2015.11.025




### DATA INPUT INSTRUCTIONS FOR MODFLOW-SDA ###
Input to the Stream Depletion Analysis (SDA) Package is read from the file that has type 'SDA' in the name file.
All variables are free format.

```
Dataset 1:  NSEN NCBMAX
Dataset 2:  FileCof
Dataset 3a: ScenName
Dataset 3b: ListFile
Dataset 3c: SFtype SFname
Dataset 3d: ZONFile
Dataset 3e: IHED, IWEL, IRCH
Dataset 3f: [HedFile]
Dataset 3g: [WelFile]
Dataset 3h: [RCHfile]
```

##### Dataset 1 #####

1. `NSEN` --- An integer value that can be specified to be positive or negative.
The absolute value of `NSEN` is equal to the total number of scenarios.
If `NSEN < 0`, SDA runs the baseline and writes the coefficients to disk.
If `NSEN > 0`, SDA reads the coefficients from disk and start the scenario analysis.
If `NSEN = 0`, SDA is inactivated.

2. `NCBMAX` --- An integer value specifying the maximum number of head-dependent boundary condition cells.

##### Dataset 2 #####

1. `FileCof` --- HDF5 File name to store baseline flow coefficients.


##### Dataset 3 #####

Repeat `NSEN` times of Dataset 3 as each set for one scenario

1. `ScenName` --- A label for the scenario run.

2. `ListFile` --- List file name for printing the running information.

3. `SFtype` --- Solver type to be used, e.g. ``PCG``, ``SIP``, ``DE4``, ``PCGN`` or ``GM``.

4. `SFname` --- File name for the solver input; the format of the solver file is identical with MODFLOW specifications.

5. `ZONFile` --- Zone file name for flow response calculation. The zone file format follows the Zone Budget specification.

6. `IHED` --- An integer value used as a flag for writing head changes.
If `IHED > 0`, head changes will be written to a binary file.
If `IHED <= 0`, head changes will not be written to disk.

7. `IWEL` --- An integer value used as a flag for well simulation in the scenario run.
If `IWEL > 0`, the well package is activated in the scenario run.
A **WEL** file will be needed. The format of the well file follows MODFLOW specification.
Please note that pumping rate used for the well package is the accretion flow compared to the baseline condition.

8. `IRCH` --- An integer value used as a flag for recharge simulation in the scenario run.
If `IRCH > 0`, the recharge package is activated in the scenario run.
A **RCH** file will be needed. The format of the recharge file follows MODFLOW specification.
Please note that recharge rate used for the recharge package is the accretion flow compared to the baseline condition.

9. `HedFile` --- File name for writing head changes. It only reads when `IHED > 0`.

10. `WelFile` --- Input file name for the well package. It only reads when `IWEL > 0`.

11. `RCHfile` --- Input file name for the recharge package. It only reads when `IRCH > 0`.


## License
SDA is distributed under the GNU Public License Version 3. For details visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
