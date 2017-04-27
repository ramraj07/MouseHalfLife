# MouseHalfLife
Single function Matlab file which can calculate half-life of injected drugs based on pharmacokinetics data. Script can handle two popular models for pharmacokinetics, mono-exponential or biexponential models.

This function is effectively a script which is used to fit pharmacokinetic
data to various exponential decay-models in order to find the half-life of
various biomolecules in mice. Input is the data from all the counts (or
ELISA values) and must be specified as a CSV (Comma Seperated Value) file
in the prescribed format. Output is a PDF with all the fitting data and
calculated half-lives.

`dataFilename= 'PkainKOandGKO03.27.10.csv';`
Enter the name of the csv file that contains the half-life data in the
prescribed format. CSV is a special data-only format you can save in excel; the
file should contain only a single table of the form:

     Group:    loxp+/+     loxp+/+     loxp+/+  ...
     Mouse:    138M-       138Mr       138Ml       ...
     0         75.4        78.6        94.3        ...
     8         64.8        66          83.4        ...
     10        61.6        63          78.6        ...
     14        58.4        58.2        71.5        ...
     24        49.3        48.7        58.6        ...
     34        42.6        39.8        45.5        ...
     54        30          27.2        30.3        ...
     71        24.4        20.6        24.4        ...
     89        18.68       15.18       17.27       ...
     110       13.44       11.08       11.76       ...
     ...       ...         ...         ...         ...
     ...       ...         ...         ...         ...
     ...       ...         ...         ...         ...

 The first column contains the timepoints at which readings were taken and
 the first two rows tell the group and mouse identification information.
 Give exactly same text for mice which are in the same group, thats how
 they are "grouped" for the report.


     ExperimentTitle = 'Pka in KO and GKO 03.27.10';
     ExperimentDescription = '';
 Provide description of the experiment, be as elaborate as you want

     ExperimentDate = '03-27-10';
     ExperimenterName=  '';

     CorrectBackground = 'no';
 If there is background correction to be done, enter 'no'. Otherwise enter
 'yes'.

     IsotopeHalfLife =1430.4;
In hours (or whatever time unit used in the data

     TimeUnits = 'hours';

     AutomaticInitialConditions = 'yes';
 If you want initial conditions to be determined automatically, put 'yes'.
 Otherwise put 'no' and provide values for the conditions below. If
 single-exponential fitting method is used only the HalfLife1 and Coeff1
 values will be utilized.

     InitialConditions_HalfLife1 = 20;
     InitialConditions_HalfLife2 = 220;
     InitialConditions_Coeff1 = 0.03;
     InitialConditions_Coeff2 = 0.003;


 BEST DEFAULT FITTING METHOD: 2
 
     fittingMethod = 1;
       
 There are three fitting methods you can use, and you choose by putting
 1,2 or 3 above. The explanations for them are:

 1 = Independent-coefficient biexponential fit 
 
       In this method the equation used is,
           y = C1*exp(-k1*t)+C2*exp(-k2*t)
       This is the traditionally used method, where you expect the curve
       to follow a biexponential curve. k1 and k2 are the half lives.

 2 = Semi-dependent-coefficient biexponential fit
 
       The equation used in this method is,
           y = 0.001*(C1*C2*exp(-k1*t)+C1*(100-C2)exp(-k2*t))
       While this equation might look different, it means exactly the same
       as the Independent-coefficient equation; the only difference is
       that here the Constants C1 and C2 indicate physiologically relevant
       values, viz. C1 is the approximate starting value at which the
       readings start and C2 is the  contribution of the first
       exponential to the whole decay process. The results obtained from
       (1) or (2) are generally similar upto the third digit (after that
       the decimals change a bit because of nuances in the fitting
       landscape)

 3 = Single exponential fit
 
       The equation used here is,
           y = C*exp(-k*t)
       This should be used in special cases where the decay looks more
       like a single exponential one rather than biexponential
