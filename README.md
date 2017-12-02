# mf2005-sda
## A MODFLOW package to linearize Stream Depletion Analysis ##

### DATA INPUT INSTRUCTIONS FOR MODFLOW-SDA ###
Input to the Stream Depletion Analysis (SDA) Package is read from the file that has type 'SDA' in the name file. 
All variables are free format.

```
Dataset 1:  NSEN NCBMAX
Dataset 2a: FileCof
Dataset 2b: FileOut
Dataset 3a: ScenName
Dataset 3b: ListFile
Dataset 3c: SFtype SFname
Dataset 3d: ZONFile
Dataset 3e: IWEL, IRCH, IHED
Dataset 3f: [WelFile]
Dataset 3g: [RCHfile]
Dataset 3h: [HedFile]
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

2. `FileOut` --- HDF5 File name to store results of scenario runs. 
If this is a blank line, it will write to indivdual CSV2 files named by `ScenName`.

##### Dataset 3 #####

Repeat `NSEN` times of Dataset 3 as each set for one scenario

1. `ScenName` --- A label for the scenario run. 

2. `ListFile` --- List file name for printing the running information. 

3. `SFtype` --- Solver type to be used, e.g. ``PCG", ``SIP", ``DE4", ``PCGN" or ``GMG". 

4. `SFname` --- File name for the solver input; the format of the solver file is identical with MODFLOW specifications.

5. `ZONFile` --- Zone file name for flow response calculation. The zone file format follows the Zone Budget specification.

6. `IWEL` --- An integer value used as a flag for well simulation in the scenario run. 
If `IWEL > 0`, the well package is activated in the scenario run. 
A **WEL** file will be needed. The format of the well file follows MODFLOW specification. 
Please note that pumping rate used for the well package is the accretion flow compared to the baseline condition.

7. `IRCH` --- An integer value used as a flag for recharge simulation in the scenario run. 
If `IRCH > 0`, the recharge package is activated in the scenario run. 
A **RCH** file will be needed. The format of the recharge file follows MODFLOW specification. 
Please note that recharge rate used for the recharge package is the accretion flow compared to the baseline condition.

8. `IHED` --- An integer value used as a flag for writing head changes. 
If `IHED < 0`, head changes will be written to the HDF5 file. 
If `IHED > 0`, head changes will be written to a formatted file with the format specification ``(10(1X1PE13.5))". 
If `IHED = 0`, head changes will not be written to disk.

9. `WelFile` --- Input file name for the well package. It only reads when `IWEL > 0`.

10. `RCHfile` --- Input file name for the recharge package. It only reads when `IRCH > 0`.

11. `HedFile` --- File name for writing head changes. It only reads when `IHED > 0`.
